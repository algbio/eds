#ifndef PBWT_HPP
#define PBWT_HPP

#include <vector>
#include <iostream>
#include "segment.hpp"
#include "msa_chunker.hpp"

using std::vector;
using std::cerr, std::endl;

// #define PBWT_DEBUG

/*
 * data structure for computing meaningful left extensions in linear time
 */
namespace pbwt {
  typedef msa_chunker::msa_chunker msa_t;
  typedef segment::seg_index seg_index;

  class pbwt {
    private:
      pbwt(const seg_index r) {
        #ifdef PBWT_DEBUG
          cerr << "Init pbwt arrays of size " << r << endl;
        #endif
      }
      ~pbwt() = default;
      // Arrays for current column
      vector<unsigned long long> ak;
      vector<unsigned long long> sk;
      vector<unsigned long long> tk;
      vector<unsigned long long> ek;
      // Extra space
      vector<unsigned long long> cnt;
      vector<unsigned long long> prev;
      vector<unsigned long long> a;
      vector<unsigned long long> e;
    public:
      // Singleton pattern
      static pbwt& instance(const seg_index r) {
        static pbwt inst(r);
        return inst;
      }
      // Delete copy and move operations
      pbwt(const pbwt&) = delete;
      pbwt& operator=(const pbwt&) = delete;
      pbwt(pbwt&&) = delete;
      pbwt& operator=(pbwt&&) = delete;
      // Column streaming algorithm
      vector<pair<seg_index, seg_index>> compute_meaningful_extensions(
            msa_t &idx,
            const seg_index r,
            const seg_index c,
            const seg_index L,
            const seg_index U,
            const seg_index y
          ) { 
        vector<pair<seg_index, seg_index>> L_y;  // 1-based indexing
        if (y < L) {
          return L_y; // No extension possible
        }
        return L_y;
      }
  };
}

#endif