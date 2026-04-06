/*
 * load a FASTA MSA with htslib into memory chunk by chunk, but the chunking
 *   is done automatically, hoping for a linear forward scan of the MSA
 */
#ifndef MSA_CHUNKER_H
#define MSA_CHUNKER_H

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>
#include <filesystem>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>

using std::vector;
using std::string;
using std::cerr, std::endl, std::flush;
using std::min, std::max;
using std::filesystem::path, std::filesystem::exists, std::filesystem::last_write_time, std::filesystem::remove;

namespace msa_chunker {
  constexpr int MAX_CHUNK_COLS = 131072;
  typedef hts_pos_t msa_pos_t;
  class msa_chunker {
    private:
      faidx_t *idx = NULL;
      msa_pos_t rows, cols, max_query_length;
      msa_pos_t chunk_start = std::numeric_limits<msa_pos_t>::max(), chunk_cols = -1;
      vector<string> msa_chunk;

      void load_chunk(msa_pos_t startcol, msa_pos_t length) {
        assert(length <= MAX_CHUNK_COLS);
        if (chunk_start <= startcol and startcol + length <= chunk_start + chunk_cols) {
          //cerr << "DEBUG: not loading chunk (query was [" << startcol << ".." << startcol + length - 1 << "])" << endl;
          return;
        }

        msa_chunk.clear();
        chunk_start = max((msa_pos_t)0, startcol - max_query_length);
        chunk_cols = min((msa_pos_t)MAX_CHUNK_COLS, cols - chunk_start);
        //cerr << "DEBUG: loading chunk [" << chunk_start << ".." << chunk_start + chunk_cols - 1 << "] (0-based) (query was [" << startcol << ".." << startcol + length - 1 << "])" << endl;

        for (msa_pos_t r = 0; r < rows; r++) {
          msa_pos_t out_len;
          char *s_str = faidx_fetch_seq64(idx, faidx_iseq(idx, r), chunk_start, chunk_start + chunk_cols - 1, &out_len);
          assert(out_len == chunk_cols);
          msa_chunk.push_back(string(s_str));
          free(s_str);
        }
      }

    public:
      msa_chunker() = delete;
      msa_chunker(const string &fastapath, const int max_qlen) {
        max_query_length = max_qlen;
        path fastap(fastapath);
        path fastaindex(fastapath + ".fai");
        path gzip_index(fastapath + ".gzi");
        bool compressed = false;
        {
          BGZF *bgzf = bgzf_open(fastap.c_str(), "r");
          if (bgzf->is_compressed) {
            compressed = true;
          }
          bgzf_close(bgzf);
        }

        if (!exists(fastaindex) or (compressed and !exists(gzip_index))) {
          cerr << "Index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\"" : "") << " not found, generating the index" << ((compressed) ? "es" : "") << "..." << flush;
          if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
            cerr << "\nERROR: failed to create index" << endl;
            exit(1);
          }
          cerr << " done." << endl;
        } else if (last_write_time(fastaindex) < last_write_time(fastap) or
            (compressed and last_write_time(gzip_index) < last_write_time(fastap))) {
          cerr << "Index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\" are" : " is") << " older than MSA, regenerating the index" << ((compressed) ? "es" : "") << "..." << flush;
          remove(fastaindex);
          if (compressed and exists(gzip_index))
            remove(gzip_index);
          if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
            cerr << "\nERROR: failed to create index" << endl;
            exit(1);
          }
          cerr << " done." << endl;
        } else {
          cerr << "Index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " and \"" + gzip_index.string() + "\"" : "") << " found" << endl;
        }

        assert(exists(fastaindex) and (!compressed or exists(gzip_index)));
        if (!(idx = fai_load3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str(), FAI_NONE))) {
          cerr << "\nERROR: failed to create index" << endl;
          exit(1);
        }

        rows = faidx_nseq(idx);
        int c = -1;
        for (int i = 0; i < rows; i++) {
          const char* seq_name = faidx_iseq(idx, i);
          const msa_pos_t seq_len = faidx_seq_len64(idx, seq_name);
          if (c == -1) {
            c = seq_len;
          } else if (seq_len != c) {
            cerr << "ERROR: MSA has rows of different length! (" << string(seq_name) << ")" << endl;
            exit(1);
          }
        }
        cols = c;
      }

      msa_pos_t get_rows() {
        return rows;
      }
      msa_pos_t get_cols() {
        return this->cols;
      }

      string msa_substr(msa_pos_t row, msa_pos_t col, msa_pos_t length) {
        assert(0 <= row and row < rows and 0 <= col and col < cols and col + length <= cols);
        load_chunk(col, length);
        return msa_chunk[row].substr(col - chunk_start, length);
      }

      ~msa_chunker() {
        fai_destroy(idx);
      }
  };
} // msa_chunker
#endif
