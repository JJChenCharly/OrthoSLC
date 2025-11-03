# OrthoSLC Developer Guide

This guide explains how the orchestration scripts and helper libraries in `src/` fit together.  It
is intended to let you understand the control flow, major data structures, and concurrency model
without reading every file line by line.

## 1. High-level workflow

The core workflow is a nine-step pipeline that starts from raw genome files and finishes with
clustered gene families.  Each numbered step has a dedicated executable under `src/` that expects
tabular or FASTA outputs from the previous step.  The diagram below captures the dependencies:

```
Step1_preWalk.py           → strain file manifest (TSV)
Step2_simple_derep.cpp     → dereplicated FASTA files (+ copy map)
Step3_pre_cluster.cpp      → concatenated derep FASTA + metadata tables + species bins
Step4_makeblastdb.py       → BLAST databases per strain
Step4_reciprocal_blast.py  → pairwise BLAST tables
Step5_query_binning.cpp    → hashed bins of BLAST hits
Step6_filter_n_bin.cpp     → filtered, weighted edges per bin
Step7_RBF.cpp              → reciprocal best-friend (RBF) edges
Step8_MCL.py or Step8_SLC.cpp → graph clustering (MCL or streaming)
cluster_fusion.cpp         → merge MCL clusters with pre-clusters
Step9_write_clusters.cpp   → final FASTA/TSV outputs per cluster type
```

Supporting modules such as `OrthoSLC.py`, `ThreadPool.h`, `Utils.hpp/.cpp`, and the toolkit scripts
(`TK_*.py`) provide shared infrastructure and optional downstream analyses.

## 2. Step-by-step behaviour

### Step 1 – Collect inputs (`Step1_preWalk.py`)
* Recursively walks an input directory, collecting every file that matches the requested file
  extension and writing a manifest that records a numeric index, a “parent” identifier inferred from
the containing folder, and the absolute path.【F:src/Step1_preWalk.py†L6-L72】
* Validates that the input directory and the output TSV’s parent directory exist before writing the
manifest.【F:src/Step1_preWalk.py†L52-L66】

### Step 2 – Dereplication (`Step2_simple_derep.cpp`)
* Loads the Step 1 TSV into a map, then uses a thread pool to deduplicate each strain’s FASTA file.
  Unique sequences are written to strain-specific FASTA files while duplicate IDs are recorded in a
  shared copy-info map protected by a mutex.【F:src/Step2_simple_derep.cpp†L13-L105】
* Emits two artefacts: dereplicated FASTA files (one per strain) and a TSV that traces duplicate
  accession IDs back to the representative sequence.【F:src/Step2_simple_derep.cpp†L107-L142】

### Step 3 – Pre-clustering (`Step3_pre_cluster.cpp`)
* Concatenates all dereplicated FASTA records, computes per-sequence length and ID metadata, and
  groups sequences by species prefix (text before the first dash).【F:src/Step3_pre_cluster.cpp†L13-L107】
* Writes three products: a combined FASTA for downstream BLAST, per-strain “nr” FASTA files, and
  tabular summaries (`seq_len_info`, `id_info`, and `pre_cluster`). Parallel writing uses
  `ThreadPool` to accelerate file emission.【F:src/Step3_pre_cluster.cpp†L108-L144】
* Stores only pointers to the deduplicated `std::string` instances returned by `DeduplicatorExact`
  when building the per-species maps, so each FASTA sequence is held once in memory and reused by
  the writers that materialise species-level NR files.【F:src/Step3_pre_cluster.cpp†L115-L147】

### Step 4a – Build BLAST databases (`Step4_makeblastdb.py`)
* Wraps the `makeblastdb` executable in the `makeblastdb` tasker from `OrthoSLC.py`, creating an
  output directory with one database per dereplicated FASTA.【F:src/Step4_makeblastdb.py†L5-L73】
* Schedules database creation through the pseudo-pool interface defined in `OrthoSLC.py`, using the
  requested number of worker processes.【F:src/Step4_makeblastdb.py†L75-L107】

### Step 4b – Reciprocal BLAST (`Step4_reciprocal_blast.py`)
* Measures each BLAST database’s size, sorts jobs by descending size, and feeds them to the dynamic
  scheduler so larger databases receive more CPU threads.【F:src/Step4_reciprocal_blast.py†L63-L126】
* Executes BLAST (`blastn` or `blastp`) with user-specified parameters, writing tabular results into
  an output directory mirroring the database structure.【F:src/Step4_reciprocal_blast.py†L128-L175】

### Step 5 – Query binning (`Step5_query_binning.cpp`)
* Hashes each BLAST hit into one of `bin_level` buckets based on the interacting pair, enabling
  downstream processing to operate on smaller files. Each BLAST output file is processed in parallel
  with thread-safe writers (optional `--no_lock_mode` bypasses mutexes when using isolated files).【F:src/Step5_query_binning.cpp†L1-L117】【F:src/Step5_query_binning.cpp†L118-L183】
* Tune `--bin_level` to balance resource usage: higher values spread hits across more files so each
  `bins_to_save[b]` vector stays smaller in memory, while lower values reduce file-count overhead but
  produce larger per-bin buffers and outputs. Pick a level that matches available RAM and I/O
  bandwidth for the BLAST workload.【F:src/Step5_query_binning.cpp†L13-L86】【F:src/Step5_query_binning.cpp†L128-L175】

### Step 6 – Filter and weight bins (`Step6_filter_n_bin.cpp`)
* Reads the sequence length table and pre-cluster assignments to filter BLAST bins by length ratio,
  similarity threshold, and desired weight field (`bitscore`, `pident`, or `evalue`).【F:src/Step6_filter_n_bin.cpp†L1-L91】
* Emits weighted edges into per-bin files, again using mutex-protected writers unless `--no_lock`
  is enabled.【F:src/Step6_filter_n_bin.cpp†L92-L173】

### Step 7 – Reciprocal best friends (`Step7_RBF.cpp`)
* Consumes the filtered bins, keeping only reciprocal best-hit pairs per species combination. The
  job distribution mirrors Step 5/6, with optional lock-free mode when safe.【F:src/Step7_RBF.cpp†L1-L96】【F:src/Step7_RBF.cpp†L97-L173】

### Step 8 – Graph clustering
* **MCL route (`Step8_MCL.py`):** Streams all Step 7 edge files into an external `mcl` command and
  pipes the clustering back into `cluster_fusion`, which merges the results with Step 3
  pre-clusters before writing `0.txt` in the target directory.【F:src/Step8_MCL.py†L1-L79】【F:src/Step8_MCL.py†L81-L116】
* **Streaming large-cluster route (`Step8_SLC.cpp`):** Batches Step 7 outputs into compression
  groups, optionally appending the Step 3 pre-cluster file, and iteratively merges clusters without
  invoking `mcl`. Supports an “all-in-one” mode (`--compression_size all`) that processes every
  input and the pre-cluster file together.【F:src/Step8_SLC.cpp†L1-L120】【F:src/Step8_SLC.cpp†L121-L207】
* `cluster_fusion.cpp` provides the merger used by the MCL route: it loads the Step 3 pre-cluster
  sets, unions them with each MCL cluster received on stdin, and ensures orphan genes still appear in
  the final output.【F:src/cluster_fusion.cpp†L1-L80】

### Step 9 – Materialise cluster outputs (`Step9_write_clusters.cpp`)
* Loads final clusters, dereplicated FASTA, and metadata tables in parallel, then writes cluster
  FASTA files and summary TSVs per requested cluster type (accessory, strict core, surplus core).
  Optional percentage thresholds restrict accessory clusters to those shared across a minimum
  fraction of genomes.【F:src/Step9_write_clusters.cpp†L1-L108】【F:src/Step9_write_clusters.cpp†L109-L173】

## 3. Shared infrastructure

### `OrthoSLC.py`
* Defines reusable “tasker” classes that know how to build command lists for `makeblastdb` and BLAST
  runs while sharing scheduling logic. The `Psudo_pool_excute` method dynamically assigns more
  threads to larger jobs, whereas `Pool_excute` wraps Python’s `multiprocessing.Pool` for simpler
  workloads.【F:src/OrthoSLC.py†L1-L157】【F:src/OrthoSLC.py†L158-L233】

### `ThreadPool.h`
* Provides a lightweight C++ thread pool with futures support. It is used throughout the C++ steps
  to distribute file-level work items while keeping writers thread-safe via external mutexes.【F:src/ThreadPool.h†L1-L83】

### `Utils.hpp` / `Utils.cpp`
* Implements streaming FASTA readers/writers and deduplication helpers (`DeduplicatorHash` and
  `DeduplicatorExact`), along with small TSV utilities for reading manifests, logging sequence
  lengths, and writing metadata tables.【F:src/Utils.hpp†L1-L125】【F:src/Utils.cpp†L1-L79】

## 4. Alignment and SNP toolkit scripts

These optional scripts operate on cluster outputs or alignments produced outside the main pipeline.

* **`TK_mafft.py` / `TK_kalign.py`:** Batch wrappers around the MAFFT and Kalign aligners. They reuse
  the pseudo-pool scheduler to prioritise larger input FASTA files and stream aligned outputs to an
  output directory, optionally passing extra command-line arguments straight through to the aligner.【F:src/TK_mafft.py†L1-L119】
* **`TK_AlnConcat.py`:** Concatenates per-gene alignments into a supermatrix, using the Step 1 ID map
  to ensure records stay in a consistent strain order before writing in the requested format.【F:src/TK_AlnConcat.py†L1-L66】
* **`TK_SNPmat.py`:** Builds a pairwise SNP count matrix across all alignments by converting sequences
  to NumPy byte arrays, computing upper-triangular differences in parallel, and remapping results
  back to strain names via the Step 1 ID TSV.【F:src/TK_SNPmat.py†L1-L92】【F:src/TK_SNPmat.py†L93-L151】

## 5. Data hand-off checklist

When wiring new functionality into the pipeline, ensure the following contracts hold:

* Filenames produced by Step 2 (`<strain>.fasta`) match the strain identifiers recorded in the Step 1
  TSV, because later stages derive species names by splitting on the first dash.【F:src/Step2_simple_derep.cpp†L66-L105】
* Step 3’s `pre_cluster` and `id_info` tables must accompany any rerun of Steps 8–9; both the MCL
  fusion tool and the final writer rely on them for expanding clusters and recovering original
  headers.【F:src/Step3_pre_cluster.cpp†L108-L144】【F:src/Step9_write_clusters.cpp†L1-L108】
* Whenever you adjust binning parameters (Steps 5–7), update downstream scripts to consume the same
  bin level and lock-mode expectations to avoid mismatched file counts.【F:src/Step5_query_binning.cpp†L1-L117】【F:src/Step7_RBF.cpp†L1-L96】

By following the relationships documented above, you can extend or troubleshoot each stage of the
OrthoSLC workflow without diving into the underlying implementation details.
