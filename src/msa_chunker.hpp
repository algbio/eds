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
    protected:
      constexpr static msa_pos_t MIN_CHUNK_COLS = 131072;
      msa_pos_t rows, cols, max_chunk_cols;
      msa_pos_t chunk_start = std::numeric_limits<msa_pos_t>::max(), chunk_cols = -1;
      vector<string> msa_rows;
      vector<string> msa_cols;
      // Which chunks to store in memory
      bool row_major = true;
      bool column_major = true;

      virtual void load_chunk(msa_pos_t startcol, msa_pos_t length) = 0;
      virtual string msa_substr_long(msa_pos_t row, msa_pos_t col, msa_pos_t length) = 0;
    public:
      virtual ~msa_chunker() = default;
      msa_pos_t get_rows() {
        return rows;
      }
      msa_pos_t get_cols() {
        return cols;
      }
      void set_row_major(bool rm) {
        row_major = rm;
        msa_rows.clear();
        chunk_start = std::numeric_limits<msa_pos_t>::max();
        chunk_cols = -1;
      }
      void set_column_major(bool cm) {
        column_major = cm;
        msa_cols.clear();
        chunk_start = std::numeric_limits<msa_pos_t>::max();
        chunk_cols = -1;
      }
      // Return the sequence at MSA[row, col..col+length]
      string msa_substr(msa_pos_t row, msa_pos_t col, msa_pos_t length) {
        assert(0 <= row and row < rows and 0 <= col and col < cols and col + length <= cols);
        assert (row_major);
        if (length < max_chunk_cols) {
          load_chunk(col, length);
          return msa_rows[row].substr(col - chunk_start, length);
        } else {
          // Not in chunk memory
          return msa_substr_long(row, col, length);
        }
      }
      // Return the character at MSA[row, col]
      char msa_at(msa_pos_t row, msa_pos_t col) {
        assert(0 <= row and row < rows and 0 <= col and col < cols);
        assert(column_major or row_major);
        load_chunk(col, 1);
        if(column_major){
          return msa_cols[col - chunk_start][row];
        } else {
          return msa_rows[row][col - chunk_start];
        }
      }
      // Return a column from the MSA
      string& get_column(msa_pos_t col) {
        assert (0 <= col and col < cols);
        assert (column_major);
        load_chunk(col, 1);
        return msa_cols[col - chunk_start];
      }
  };

  /*
   * load a FASTA MSA with htslib into memory chunk by chunk, but the chunking
   *   is done automatically, hoping for a linear forward scan of the MSA; since
   *   we expect all rows to be processed, we only consider vertical MSA slices
   */
  class fasta_chunker : public msa_chunker {
    private:
      faidx_t *idx = NULL;

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

        msa_rows.clear();
        chunk_start = startcol;
        chunk_cols = min(max_chunk_cols, cols - chunk_start);
#ifdef MSA_CHUNKER_DEBUG
        cerr << "DEBUG: loading chunk [" << chunk_start << ".." << chunk_start + chunk_cols - 1 << "] (0-based) (query was [" << startcol << ".." << startcol + length - 1 << "])" << endl;
#endif

        // Calculate row-major chunk
        if (row_major) {
          for (msa_pos_t r = 0; r < rows; ++r) {
            msa_pos_t out_len;
            char *s_str = faidx_fetch_seq64(idx, faidx_iseq(idx, r), chunk_start, chunk_start + chunk_cols - 1, &out_len);
            assert(out_len == chunk_cols);
            msa_rows.push_back(string(s_str));
            free(s_str);
          }
        }
      
        // Calculate column-major chunk
        if (column_major) {
          // Relies on msa_rows
          assert(row_major);
          msa_cols.resize(chunk_cols);
          for (msa_pos_t j = 0; j < chunk_cols; j++) {
            string column = "";
            for (msa_pos_t i = 0; i < rows; i++) {
              column.push_back(msa_rows[i][j]);
            }
            msa_cols[j] = column;
          }
        }
      }

      string msa_substr_long(msa_pos_t row, msa_pos_t col, msa_pos_t length) override {
        msa_pos_t out_len;
        char *s_str = faidx_fetch_seq64(idx, faidx_iseq(idx, row), col, col + length - 1, &out_len);
        assert(out_len == length);
        string s(s_str);
        free(s_str);
        return s;
      }

    public:
      fasta_chunker() = delete;

      /*
       * index a given (gzipped) FASTA file
       */
      fasta_chunker(const string &fastapath, const msa_pos_t max_qlen, const bool verbose) {
        max_chunk_cols = max(max_qlen, MIN_CHUNK_COLS);
        if(!std::filesystem::is_regular_file(fastapath)){
          throw runtime_error("ERROR: FASTA file could not be found");
        }
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

      ~fasta_chunker() {
        fai_destroy(idx);
      }
  };

  /*
   * load a column major msa into memory chunk by chunk, should speed up column streaming
   */
  class column_chunker : public msa_chunker {
    private:
      std::streampos matrix_start;
      ifstream msa_file;
      
      /*
       * explicitly load chunk [startcol..startcol+length) into memory
       */
      void load_chunk(msa_pos_t startcol, msa_pos_t length) {
        assert(length <= max_chunk_cols);
        if (chunk_start <= startcol and startcol + length <= chunk_start + chunk_cols) {
          return;
        }
        
        // Read chunk from transposed msa file
        chunk_start = startcol;
        chunk_cols = min(max_chunk_cols, cols - chunk_start);
        std::streampos col_pos = matrix_start + chunk_start * (rows + 1);
        msa_file.seekg(col_pos);
        vector<char> msa_buffer(chunk_cols * (rows + 1));
        msa_file.read(msa_buffer.data(), chunk_cols * (rows + 1));

        // Split the column-major chunk
        if(column_major) {
          msa_cols.resize(chunk_cols);
          for (msa_pos_t j = 0; j < chunk_cols; j++) {
            auto msa_col = msa_buffer.begin() + j * (rows + 1);
            msa_cols[j] = string(msa_col, msa_col + rows);
          }
        }

        // Compute original row-major chunk
        if(row_major) {
          msa_rows.resize(rows);
          for (msa_pos_t i = 0; i < rows; i++) {
            string row = "";
            for (msa_pos_t j = 0; j < chunk_cols; j++) {
              row.push_back(msa_buffer[j * (rows + 1) + i]);
            }    
            msa_rows[i] = row;
          }
        }
      }

       // Not implemented
      string msa_substr_long(msa_pos_t row, msa_pos_t col, msa_pos_t length) override {
        throw runtime_error("ERROR: substring length exceeds maximum column chunk value");
      }

    public:
      column_chunker() = delete;

      column_chunker(const string &msapath, const msa_pos_t max_qlen, const bool verbose) {
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

      ~column_chunker(){
        msa_file.close();
      }
  };
} // msa_chunker
#endif // MSA_CHUNKER_HPP
