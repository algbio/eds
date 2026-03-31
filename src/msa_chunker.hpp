/*
 * load a FASTA MSA with htslib into memory chunk by chunk, but the chunking
 *   is done automatically, hoping for a linear forward scan of the MSA
 */
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#define MAX_CHUNK_COLS 131072

using std::vector;
using std::string;
using std::cerr, std::endl;
using std::min, std::max;

namespace msa_chunker {
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
        if (!(idx = fai_load3_format(fastapath.c_str(), NULL, NULL, FAI_CREATE, FAI_FASTA))) {
          cerr << "Failed to load/create index" << endl;
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
            cerr << "MSA has rows of different length! (" << string(seq_name) << ")" << endl;
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
