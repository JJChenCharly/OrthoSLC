#include "Utils.hpp"
#include <algorithm>

namespace fasta {

Reader::Reader(const std::string& filename) : file_(filename)
{
    if (!file_) throw std::runtime_error("Unable to open FASTA file: " + filename);
}

bool Reader::next(std::string& header, std::string& seq)
{
    header.clear();
    seq.clear();

    std::string line;

    // If we have a cached header from the previous read-ahead, use it first.
    if (!cached_header_.empty()) {
        header.swap(cached_header_);
    } else {
        while (std::getline(file_, line)) {
            if (!line.empty() && line[0] == '>') {
                header = line.substr(1);
                break;
            }
        }
        if (header.empty()) return false; // EOF
    }

    std::ostringstream os;
    while (std::getline(file_, line)) {
        if (!line.empty() && line[0] == '>') {
            cached_header_ = line.substr(1); // save for the next call
            break;
        }
        os << line;
    }

    if (file_.eof()) eof_ = true;
    seq = os.str();
    return true;
}

Writer::Writer(const std::string& filename, std::size_t buffer_limit)
    : file_(filename, std::ios::trunc), buffer_limit_(buffer_limit)
{
    if (!file_) throw std::runtime_error("Unable to open output FASTA: " + filename);
}

void Writer::format_sequence(const std::string& seq, std::string& out, std::size_t line_len)
{
    const std::size_t n = seq.size();
    for (std::size_t i = 0; i < n; i += line_len) {
        out.append(seq, i, std::min(line_len, n - i));
        out.push_back('\n');
    }
}

void Writer::write(const std::string& full_header, const std::string& seq)
{
    std::string record;
    record.reserve(full_header.size() + seq.size() + seq.size() / 60 + 4);
    record.push_back('>');
    
    record.append(full_header);
    record.push_back('\n');
    format_sequence(seq, record);

    buffer_.append(record);
    if (buffer_.size() >= buffer_limit_) flush();
}

void Writer::flush()
{
    file_ << buffer_;
    buffer_.clear();
}

Writer::~Writer() { flush(); }

bool DeduplicatorHash::add(const std::string& seq, std::string id) {
    const std::size_t h = hasher_(seq);
    // Try to insert an empty vector if 'h' is not present
    auto it = map_.try_emplace(h, IdList{}).first;
    // Move 'id' into the vector to avoid an extra copy
    it->second.emplace_back(std::move(id));
    return it->second.size() == 1;  // true if first occurrence
}
std::pair<const std::size_t*, bool> DeduplicatorHash::add_and_get(const std::string& seq, const std::string& id) {
    const std::size_t h = hasher_(seq);
    // try_emplace will move `seq` into the key on insertion, with zero extra copies
    auto [it, inserted] = map_.try_emplace(h, IdList{});
    it->second.emplace_back(std::move(id));            // record the ID
    return { &it->first, inserted };
}

bool DeduplicatorExact::add(std::string seq, std::string id) {
    // Move 'seq' into the key to avoid a copy
    auto it = map_.try_emplace(std::move(seq), IdList{}).first;
    // Move 'id' into the vector
    // because id never replicate
    it->second.emplace_back(std::move(id));
    return it->second.size() == 1;
}
std::pair<const std::string*, bool> DeduplicatorExact::add_and_get(const std::string& seq, const std::string& id) {
    // try_emplace will move `seq` into the key on insertion, with zero extra copies
    auto [it, inserted] = map_.try_emplace(std::move(seq), IdList{});
    it->second.emplace_back(std::move(id));            // record the ID
    return { &it->first, inserted };
}

} // namespace fasta

//======================  example.cpp  ======================
/*
#include "FastaIO.h"
using namespace fasta;

int main() {
    Reader reader("input.fasta");
    Writer writer("unique.fasta", 2 * 1024 * 1024);
    DeduplicatorHash dedup;  // or DeduplicatorExact for perfect identity

    std::string header, seq;
    while (reader.next(header, seq)) {
        const std::string id = header.substr(0, header.find(' '));
        if (dedup.add(seq, id)) {
            writer.write(id, header, seq);
        }
    }
    writer.flush();
}
*/
