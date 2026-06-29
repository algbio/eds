#ifndef MSA_CHUNKER_HPP
#define MSA_CHUNKER_HPP

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <chrono>
#include <iomanip>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>

using std::vector;
using std::string, std::to_string;
using std::cerr, std::endl, std::flush;
using std::min, std::max;
using std::chrono::high_resolution_clock, std::chrono::duration_cast, std::chrono::milliseconds;
using std::setprecision;
using std::filesystem::path, std::filesystem::exists, std::filesystem::last_write_time, std::filesystem::remove;
using std::runtime_error;
using std::ifstream;

//#define MSA_CHUNKER_DEBUG

namespace msa_chunker {
  typedef hts_pos_t msa_pos_t; // type for MSA cols, rows index

  /*
   * msa can be read into chunks either through a fasta file or a txt file (column major msa)
   */
  class msa_chunker {
    public:
      virtual ~msa_chunker() = default;
      virtual msa_pos_t get_rows() = 0;
      virtual msa_pos_t get_cols() = 0;
      virtual string msa_substr(msa_pos_t row, msa_pos_t col, msa_pos_t length) = 0;
      virtual char msa_at(msa_pos_t row, msa_pos_t col) = 0;
  };

  /*
   * load a FASTA MSA with htslib into memory chunk by chunk, but the chunking
   *   is done automatically, hoping for a linear forward scan of the MSA; since
   *   we expect all rows to be processed, we only consider vertical MSA slices
   */
  class fasta_chunker : public msa_chunker {
    private:
      constexpr static msa_pos_t MIN_CHUNK_COLS = 131072;

      faidx_t *idx = NULL;
      msa_pos_t rows, cols, max_chunk_cols;
      msa_pos_t chunk_start = std::numeric_limits<msa_pos_t>::max(), chunk_cols = -1;
      vector<string> msa_chunk;

      /*
       * explicitly load chunk [startcol..startcol+length) into memory
       */
      void load_chunk(msa_pos_t startcol, msa_pos_t length) {
        assert(length <= max_chunk_cols);
        if (chunk_start <= startcol and startcol + length <= chunk_start + chunk_cols) {
#ifdef MSA_CHUNKER_DEBUG
          cerr << "DEBUG: not loading chunk (query was [" << startcol << ".." << startcol + length - 1 << "])" << endl;
#endif
          return;
        }

        msa_chunk.clear();
        chunk_start = startcol;
        chunk_cols = min(max_chunk_cols, cols - chunk_start);
#ifdef MSA_CHUNKER_DEBUG
        cerr << "DEBUG: loading chunk [" << chunk_start << ".." << chunk_start + chunk_cols - 1 << "] (0-based) (query was [" << startcol << ".." << startcol + length - 1 << "])" << endl;
#endif

        for (msa_pos_t r = 0; r < rows; ++r) {
          msa_pos_t out_len;
          char *s_str = faidx_fetch_seq64(idx, faidx_iseq(idx, r), chunk_start, chunk_start + chunk_cols - 1, &out_len);
          assert(out_len == chunk_cols);
          msa_chunk.push_back(string(s_str));
          free(s_str);
        }
      }

    public:
      fasta_chunker() = delete;

      /*
       * index a given (gzipped) FASTA file
       */
      fasta_chunker(const string &fastapath, const msa_pos_t max_qlen, const bool verbose) {
        max_chunk_cols = max(max_qlen, MIN_CHUNK_COLS);
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
          auto start = high_resolution_clock::now();
          cerr << "Index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\"" : "") << " not found, generating the index" << ((compressed) ? "es" : "") << "..." << flush;
          if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
            cerr << "\nERROR: failed to create index" << endl;
            exit(1);
          }
          auto stop = high_resolution_clock::now();
          auto duration = duration_cast<milliseconds>(stop - start);
          cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
        } else if (last_write_time(fastaindex) < last_write_time(fastap) or
            (compressed and last_write_time(gzip_index) < last_write_time(fastap))) {
          auto start = high_resolution_clock::now();
          cerr << "Index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\" are" : " is") << " older than MSA, regenerating the index" << ((compressed) ? "es" : "") << "..." << flush;
          remove(fastaindex);
          if (compressed and exists(gzip_index))
            remove(gzip_index);
          if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
            cerr << "\nERROR: failed to create index" << endl;
            exit(1);
          }
          auto stop = high_resolution_clock::now();
          auto duration = duration_cast<milliseconds>(stop - start);
          cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
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

      msa_pos_t get_rows() override {
        return rows;
      }

      msa_pos_t get_cols() override {
        return cols;
      }

      string msa_substr(msa_pos_t row, msa_pos_t col, msa_pos_t length) override {
        assert(0 <= row and row < rows and 0 <= col and col < cols and col + length <= cols);
        if (length <= max_chunk_cols) {
          load_chunk(col, length);
          return msa_chunk[row].substr(col - chunk_start, length);
        } else {
          msa_chunk.clear();
          chunk_start = std::numeric_limits<msa_pos_t>::max();
          chunk_cols = -1;

          msa_pos_t out_len;
          char *s_str = faidx_fetch_seq64(idx, faidx_iseq(idx, row), col, col + length - 1, &out_len);
          assert(out_len == length);
          string s(s_str);
          free(s_str);
          return s;
        }
      }

      // Return the character at MSA[row, col]
      char msa_at(msa_pos_t row, msa_pos_t col) override {
        assert(0 <= row and row < rows and 0 <= col and col < cols);
        load_chunk(col, 1);
        return msa_chunk[row][col - chunk_start];
      }

      ~fasta_chunker() {
        fai_destroy(idx);
      }
  };

  /*
   * load a column major msa into memory chunk by chunk, should speed up column streaming
   */
  class column_chunker : public msa_chunker {
    private:
      constexpr static msa_pos_t MIN_CHUNK_COLS = 131072;

      msa_pos_t rows, cols, max_chunk_cols;
      msa_pos_t chunk_start = std::numeric_limits<msa_pos_t>::max(), chunk_cols = -1;
      vector<char> msa_transposed_chunk;
      vector<char> msa_chunk; // used for faster msa_substr
      std::streampos matrix_start;
      ifstream msa_file;
      
      /*
       * explicitly load chunk [startcol..startcol+length) into memory
       */
      void load_chunk(msa_pos_t startcol, msa_pos_t length, bool transpose) {
        assert(length <= max_chunk_cols);
        if (chunk_start <= startcol and startcol + length <= chunk_start + chunk_cols) {
          return;
        }
        
        // Read chunk from transposed msa file
        chunk_start = startcol;
        chunk_cols = min(max_chunk_cols, cols - chunk_start);
        std::streampos col_pos = matrix_start + chunk_start * (rows + 1);
        msa_file.seekg(col_pos);
        msa_transposed_chunk.resize(chunk_cols * (rows + 1));
        msa_file.read(msa_transposed_chunk.data(), chunk_cols * (rows + 1));

        // Compute original row-major chunk
        if(transpose || chunk_start == 0) {
          msa_chunk.resize(rows * chunk_cols);
          for (msa_pos_t i = 0; i < rows; i++) {
            for (msa_pos_t j = 0; j < chunk_cols; j++) {
              msa_chunk[i * chunk_cols + j] = msa_transposed_chunk[j * (rows + 1) + i];
            }    
          }
        }
      }

    public:
      column_chunker() = delete;

      column_chunker(const string &msapath, const msa_pos_t max_qlen) {
        max_chunk_cols = max(MIN_CHUNK_COLS, max_qlen);
        msa_file = ifstream(msapath, std::ios::binary);
        if (!msa_file) {
          throw runtime_error("ERROR: msa file could not be opened");
        }
        // Read msa dimensions
        string header, extra;
        getline(msa_file, header);
        std::istringstream iss(header);
        if (!(iss >> cols >> rows) || (iss >> extra)) {
          throw runtime_error("ERROR: first line must contain exactly two integers: columns and rows");
        }
        matrix_start = msa_file.tellg();
      }

      msa_pos_t get_rows() override {
        return rows;
      }
      
      msa_pos_t get_cols() override {
        return cols;
      }
      
      // Return the sequence at MSA[row, col..col+length]
      string msa_substr(msa_pos_t row, msa_pos_t col, msa_pos_t length) override {
        assert(0 <= row and row < rows and 0 <= col and col < cols and col + length <= cols);
        load_chunk(col, length, true);
        auto msa_row = msa_chunk.begin() + row * chunk_cols;
        return string(msa_row + col - chunk_start, msa_row + col - chunk_start + length);
      }

      // Return the character at MSA[row, col]
      char msa_at(msa_pos_t row, msa_pos_t col) override {
        assert(0 <= row and row < rows and 0 <= col and col < cols);
        load_chunk(col, 1, false);
        return msa_transposed_chunk[(col - chunk_start) * (rows + 1) + row];
      }

      ~column_chunker(){
        msa_file.close();
      }
  };
} // msa_chunker
#endif // MSA_CHUNKER_HPP
