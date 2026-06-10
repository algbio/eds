#ifndef PBWT_HPP
#define PBWT_HPP

#include <vector>
#include <utility>
#include <iostream>
#include <numeric>
#include <algorithm>
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

  const unsigned int alphabet_size = 6;
  inline unsigned int symbol_number(const char symbol) {
    switch (symbol) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      case '-': return 4;
      default: return 5; // N
    }
  }

  class pbwt {
    private:
      // Arrays for current column 1-based indexing
      vector<unsigned long long> ak;
      vector<unsigned long long> sk;
      vector<unsigned long long> tk;
      vector<unsigned long long> ek;
      // Counting sort arrays
      vector<unsigned long long> cnt;
      vector<unsigned long long> prev;
      // Extra space
      seg_index dy; // size of sk
      vector<unsigned long long> a;
      vector<unsigned long long> e;
      // Singleton constructor
      pbwt(const seg_index r): 
        ak(r + 1), sk(r + 2), tk(r + 2), ek(r + 1, 0), 
        cnt(alphabet_size + 1), prev(alphabet_size),
        dy(0), a(r + 1), e(r + 1)
      {
        // Initial sorted order 1, 2, 3..
        std::iota(ak.begin(), ak.end(), 0);
      }
      ~pbwt() = default;
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
        if (y < L) {
          return {}; // No extension possible
        }
        // Zero init counting sort arrays and frequency array
        cnt.assign(cnt.size(), 0);
        prev.assign(prev.size(), 0);
        tk.assign(tk.size(), 0);
        sk[++dy] = y;
        // Counting sort
        for (seg_index i = 1; i <= r; i++)
          cnt[symbol_number(idx.msa_at(i - 1, y - 1)) + 1]++;
        for (unsigned int i = 1; i < alphabet_size; i++)
          cnt[i] += cnt[i - 1];
        for (seg_index i = 1; i <= r; i++) {
          unsigned int b = symbol_number(idx.msa_at(ak[i] - 1, y - 1));
          cnt[b]++;
          a[cnt[b]] = ak[i];
          if (prev[b] == 0) 
            e[cnt[b]] = dy;
          else {
            // O(alphabet) amortized - very fast in practice
            e[cnt[b]] = *std::max_element(ek.begin() + prev[b] + 1, ek.begin() + i + 1);
          }
          prev[b] = i;
        }
        for (seg_index i = 1; i <= r; i++) {
          ak[i] = a[i];
        }
        // Calculate frequency array
        for (seg_index i = 1; i <= r; i++) {
          tk[e[i]]++;
        }
        // Shrink arrays sk and tk, array a acts as tmp
        seg_index j = 1;
        for (seg_index i = 1; i <= dy; i++) {
          if (tk[i] != 0) {
            a[i] = j;
            sk[j] = sk[i];
            tk[j] = tk[i];
            j++;
          }
        }
        dy = j - 1;
        // Fix ek array
        for (seg_index i = 1; i <= r; i++) {
          ek[i] = a[e[i]];
        }
        // Calculate meaningul extensions from pbwt arrays
        for(seg_index i = 1; i <= dy; i++){
          tk[dy - i] += tk[dy - i + 1];
        }
        vector<pair<seg_index, seg_index>> L_y; 
        for(seg_index i = 0; i < dy; i++){
          if(y - sk[dy - i] + 1 >= L and y - sk[dy - i] + 1 <= U)
            L_y.emplace_back(sk[dy - i], tk[dy - i]);
        }
        // Dummy element
        L_y.emplace_back(std::max((seg_index)0, y - U), -1);
        return L_y;
      }
  };
}

#endif