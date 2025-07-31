#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <memory>

namespace fasta {

/**
 * Streaming FASTA reader – mimics Biopython's SeqIO.parse.
 * Usage:
 *   fasta::Reader reader("in.fasta");
 *   std::string header, seq;
 *   while (reader.next(header, seq)) { / ... / }
 */
class Reader {
public:
    explicit Reader(const std::string& filename);

    /**
     * Load the next record.
     * @param header – string after '>'.
     * @param seq    – full (unwrapped) sequence.
     * @return true  if a record was read, false on EOF.
     */
    bool next(std::string& header, std::string& seq);

private:
    std::ifstream file_;
    std::string   cached_header_{};   // header read ahead for next() call
    bool          eof_{false};
};

/**
 * Buffered FASTA writer – wraps long sequences to 60 chars/line.
 * Flushes automatically when buffer ≥ buffer_limit or on destruction.
 */
class Writer {
public:
    explicit Writer(const std::string& filename, std::size_t buffer_limit = 4 * 1024 * 1024);
    void write(const std::string& full_header, const std::string& seq);
    void flush();
    ~Writer();

private:
    std::ofstream file_;
    std::string   buffer_;
    std::size_t   buffer_limit_;
    static void format_sequence(const std::string& seq, std::string& out, std::size_t line_len = 60);
};

/**
 * Deduplicator based on std::hash – stores collisions in a vector.
 */
class DeduplicatorHash {
public:
    using IdList = std::vector<std::string>;

    /**
     * @return true if the sequence (hash) is seen for the first time.
     */
    bool add(const std::string& seq, std::string id);
    std::pair<const std::size_t*, bool> add_and_get(const std::string& seq, const std::string& id);
    const std::unordered_map<std::size_t, IdList>& map() const { return map_; }
    void reserve(size_t n) {
        map_.reserve(n);
    }

private:
    std::unordered_map<std::size_t, IdList> map_;
    std::hash<std::string>                  hasher_;
};

/**
 * Exact-sequence deduplicator – heavier, but collision-free.
 */
class DeduplicatorExact {
public:
    using IdList = std::vector<std::string>;

    bool add(std::string seq, std::string id);
    std::pair<const std::string*, bool> add_and_get(const std::string& seq, const std::string& id);
    const std::unordered_map<std::string, IdList>& map() const { return map_; }
    void reserve(size_t n) {
        map_.reserve(n);
    }

private:
    std::unordered_map<std::string, IdList> map_;
};

} // namespace fasta

/**
 * Utilities for reading and writing tab-separated files with variable column counts.
 * The file format is not assumed to have a fixed number of columns per row.
 */
namespace tsv {
// A row represented as a vector of string tokens
using Row = std::vector<std::string>;

// Map type: key = first field of each row; value = all tokens of that row
using Map = std::unordered_map<std::string, Row>;

/**
 * Read a tab-separated file into a Map.
 * Each line is split on '\t'. Empty lines are skipped.
 * The first token of each non-empty line becomes the key, and all tokens are stored as the value vector.
 * Throws std::runtime_error if the file cannot be opened.
 */
inline void readTSV(const std::string& filename, Map& outMap) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open TSV file: " + filename);
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> tokens;
        std::string token;
        std::stringstream ss(line);
        while (std::getline(ss, token, '\t')) {
            tokens.emplace_back(token);
        }
        if (tokens.empty()) continue;
        const std::string key = tokens[0];
        outMap.emplace(key, std::move(tokens));
    }
}

/**
 * Write a Map to a tab-separated file.
 * Each entry in the map produces one line: all strings in the value vector joined by '\t'.
 * The order of rows corresponds to iteration order of the unordered_map.
 * Throws std::runtime_error if the file cannot be opened.
 */
// single helper, parameterized by whether to emit the key
template <bool EmitKey>
inline void writeTSVImpl(const std::string& filename,
                         const Map& inMap)
{
    std::ofstream file(filename, std::ios::trunc);
    if (!file) throw std::runtime_error("Unable to open TSV: " + filename);

    for (auto const& [key, row] : inMap) {
        if constexpr (EmitKey) {
            file << key;
            // include key and first of row
            for (size_t i = 0; i < row.size(); ++i)
                file << '\t' << row[i];
        }
        else {
            // just the row tokens, including row[0]
            for (size_t i = 0; i < row.size(); ++i) {
                if (i) file << '\t';
                file << row[i];
            }
        }
        file << '\n';
    }
}

// two zero-overhead wrappers
inline void writeTSV(const std::string& filename,
                     const Map& inMap)
{
    writeTSVImpl<false>(filename, inMap);
}

inline void writeTSVWithKey(const std::string& filename,
                            const Map& inMap)
{
    writeTSVImpl<true>(filename, inMap);
}

/**
 * This is specifically for trunc mode, and those writer who allows 
 * @brief Simple buffered TSV writer for two‐column output.
 *
 * Accumulates rows in a string buffer and flushes out to disk
 * once the buffer reaches buf_size.  On destruction, any remaining
 * data is flushed automatically.
 */
class BufferedTSVWriter {
public:
    BufferedTSVWriter(const std::string &filename,
                      std::size_t buf_size = 1 << 20)
      : ofs_(filename, std::ios::out /*trunc on first open*/),
        buf_size_(buf_size)
    {
        if (!ofs_) throw std::runtime_error("Failed to open TSV for writing: " + filename);
        buffer_.reserve(buf_size_);
    }

    ~BufferedTSVWriter() { flush(); }

    void writeEntry(const std::string &col1, const std::string &col2) {
        buffer_ += col1;
        buffer_.push_back('\t');
        buffer_ += col2;
        buffer_.push_back('\n');

        if (buffer_.size() >= buf_size_) {
            ofs_ << buffer_;
            buffer_.clear();
        }
    }

    void flush() {
        if (!buffer_.empty()) {
            ofs_ << buffer_;
            buffer_.clear();
        }
    }

private:
    std::ofstream ofs_;
    std::string   buffer_;
    std::size_t   buf_size_;
};

// —————————— helper to manage one-writer-per-file ——————————
template<typename... Args>
static BufferedTSVWriter& get_writer(const std::string &path, std::size_t buf_size) {
    static std::unordered_map<std::string, std::unique_ptr<BufferedTSVWriter>> map;
    auto it = map.find(path);
    if (it == map.end()) {
        auto w = std::make_unique<BufferedTSVWriter>(path, buf_size);
        it = map.emplace(path, std::move(w)).first;
    }
    return *it->second;
}

// —————————— modified len_info_TSVwriter ——————————
inline void len_info_TSVwriter(const std::string &file_path,
                               const std::string &key,
                               std::size_t        value,
                               std::size_t        buf_size = 1 << 20)
{
    auto &writer = get_writer(file_path, buf_size);
    writer.writeEntry(key, std::to_string(value));
}

// —————————— modified id_info_TSVwriter ——————————
inline void id_info_TSVwriter(const std::string &file_path,
                              const std::string &id,
                              const std::string &info,
                              std::size_t        buf_size = 1 << 20)
{
    auto &writer = get_writer(file_path, buf_size);
    writer.writeEntry(id, info);
}

} // namespace tsv

/*
 Example usage:

 #include "TabSeparatedIO.hpp"
 #include <iostream>

 int main() {
     // Read input.tsv into a map
     tsv::Map data;
     tsv::readTSV("input.tsv", data);

     // Iterate and print
     for (const auto& kv : data) {
         const auto& key = kv.first;
         const auto& row = kv.second;
         std::cout << "Key: " << key << " -> ";
         for (const auto& field : row) {
             std::cout << field << "\t";
         }
         std::cout << '\n';
     }

     // Write back to output.tsv
     tsv::writeTSV("output.tsv", data);
     return 0;
 }
*/
