#ifndef PBWT_HPP
#define PBWT_HPP

#include <vector>
#include <array>
#include <utility>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <optional>
#include "segment.hpp"
#include "msa_chunker.hpp"
#include "rmq.hpp"

using std::vector;
using std::array;
using std::cerr, std::endl;

// #define PBWT_DEBUG

/*
 * data structure for computing meaningful left extensions in linear time
 */
namespace pbwt {
  typedef msa_chunker::msa_chunker msa_t;
  typedef segment::seg_index seg_index;

  enum MaxRange {
    NAIVE,
    RECURSIVE,
    RMQ
  };

  // Range maximum algorithm
  inline constexpr enum MaxRange max_range = NAIVE;

  class pbwt {
    private:
      // Arrays for current column 1-based indexing
      vector<seg_index> ak;
      vector<seg_index> sk;
      vector<seg_index> tk;
      vector<seg_index> ek;
      // Counting sort arrays
      vector<seg_index> cnt;
      vector<seg_index> prev;
      // Extra space
      seg_index dy; // size of sk
      vector<seg_index> a;
      vector<seg_index> e;
      // Alphabet
      array<int, 256> alphabet = [] {
        array<int, 256> alph{};
        alph.fill(-1);

        alph['A'] = 0;
        alph['C'] = 1;
        alph['G'] = 2;
        alph['T'] = 3;
        alph['-'] = 4;

        return alph;
      }();
      unsigned int alphabet_size = 5;
      // Singleton constructor
      pbwt(const seg_index r): 
        ak(r + 1), sk(r + 2), tk(r + 2), ek(r + 1, 0), 
        dy(0), a(r + 2), e(r + 1)
      {
        cnt.resize(alphabet_size + 1);
        prev.resize(alphabet_size);
        // Initial sorted order 1, 2, 3..
        std::iota(ak.begin(), ak.end(), 0);
      }
      ~pbwt() = default;
      // Add new symbols to alphabet
      unsigned int symbol_number(const char symbol) {
        auto s = static_cast<unsigned char>(symbol);

        if (alphabet[s] == -1) {
          alphabet[s] = alphabet_size++;
          cnt.push_back(0);
          prev.push_back(0);
        }

        return alphabet[s];
      }
      // Recursive function for range max
      seg_index max_ek(seg_index j, seg_index i) {
        if (j != i) {
          ek[j] = std::max(ek[j], max_ek(ak[j], i));
          ak[j] = i + 1;
        }
        return ek[j];
      }
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
        // Symbol frequency
        cnt.assign(cnt.size(), 0);
        for (seg_index i = 1; i <= r; i++){
          auto c = idx.msa_at(i - 1, y - 1);
          auto s = symbol_number(c) + 1;
          cnt[s]++;
        }
        int symbols = 0;
        for (seg_index sym = 1; sym <= alphabet_size; sym++)
          if(cnt[sym] > 0)
            symbols++;
        if (y == 1 or symbols > 1) {
          // There is at least one new divergence => height change
          // Recompute pBWT arrays
          prev.assign(prev.size(), 0);
          tk.assign(tk.size(), 0);
          sk[++dy] = y;
          // RMQ for max ek
          std::optional<rmq<seg_index>> rm;
          if constexpr (max_range == RMQ) {
            rm.emplace(ek);
          }
          // Counting sort
          for (unsigned int i = 1; i < alphabet_size; i++)
            cnt[i] += cnt[i - 1];
          for (seg_index i = 1; i <= r; i++) {
            unsigned int b = symbol_number(idx.msa_at(ak[i] - 1, y - 1));
            cnt[b]++;
            a[cnt[b]] = ak[i];
            if constexpr (max_range == RECURSIVE) {
              ak[i] = i + 1;
            }
            if (prev[b] == 0) 
              e[cnt[b]] = dy;
            else {
              // Calculate the range maximum 
              if constexpr (max_range == NAIVE) {
                // O(alphabet) amortized - very fast in practice
                e[cnt[b]] = *std::max_element(ek.begin() + prev[b] + 1, ek.begin() + i + 1);
              } else if constexpr (max_range == RECURSIVE) {
                // O(log alphabet) - solution from the paper
                e[cnt[b]] = max_ek(prev[b] + 1, i);
              } else {
                // RMQ O(1) - best complexity but longer construction
                e[cnt[b]] = rm->query(prev[b] + 1, i);
              }
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
          // Calculate heights of extensions
          for(seg_index i = 1; i <= dy; i++){
            tk[dy - i] += tk[dy - i + 1];
          }
        }
        else {
          // Update last element of sk for first row divergence
          if (tk[dy] > 1){
            dy++;
          }
          sk[dy] = y;
          ek[1] = dy;
          tk[dy] = 1;
        }
        // Compute meaningful extensions
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