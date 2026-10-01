#pragma once

// Optional recording/replay for chemistry table file I/O (tables/chemistry_tables.h's
// open_chemistry_table()). Lets a multi-rank host make exactly one rank do the real
// disk I/O for every chemistry table file, then have every other rank replay from an
// in-memory copy instead of each independently re-opening (and re-contending for) the
// same files on a shared filesystem. The portable core itself stays host/MPI-agnostic:
// it never calls these, and when neither mode is active, open_chemistry_table() does a
// plain fopen exactly as before.
//
// Usage (by a host adapter, around a single Model::create() call):
//   TableIoCache cache;
//   chemistry_table_io_begin_record(&cache);   // on the one rank that reads from disk
//   ... Model::create(config) ...
//   chemistry_table_io_end();
//   // broadcast `cache` to the other ranks by whatever means (e.g. MPI)
//   chemistry_table_io_begin_replay(&cache);   // on every other rank
//   ... Model::create(config) ...               // served entirely from `cache`
//   chemistry_table_io_end();

#include <string>
#include <unordered_map>
#include <vector>

using TableIoCache = std::unordered_map<std::string, std::vector<char>>;

// Subsequent open_chemistry_table() calls do a real fopen as usual, but also capture
// each successfully-opened file's full contents into *cache (keyed by the filename
// argument, e.g. "reaction_network.tsv" -- not the full path).
void chemistry_table_io_begin_record(TableIoCache *cache);

// Subsequent open_chemistry_table() calls are served only from *cache (via fmemopen);
// a requested filename missing from the cache is a hard failure (NULL + stderr
// message), never a silent fall-through to disk.
void chemistry_table_io_begin_replay(const TableIoCache *cache);

// Returns to the default: open_chemistry_table() does a plain fopen, no recording.
void chemistry_table_io_end(void);
