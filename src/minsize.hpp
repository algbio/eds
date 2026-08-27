#ifndef MINSIZE_HPP
#define MINSIZE_HPP

#include "segment.hpp"
#include "msa_chunker.hpp"
#include "algo.hpp"
#include "RMaxQTree.h"
#include "rmqueue.h"

#include <numeric>

using std::vector, std::tie;

namespace minsize {
  typedef msa_chunker::msa_chunker msa_t;
  typedef segment::seg_index seg_index;

  // Which data structure to use for the sliding window minimum range queries
  enum MinRange {
    RING, // Naive ring buffer
    TREE, // RMaxQTree
    QUEUE // RMQueue
  };

  inline constexpr enum MinRange min_range = RING;

  // Class for solving the min-L-size segmentation problem 
  class minsize {
      msa_t &idx;
      const seg_index r;
      const seg_index c; 
      const seg_index L; 
      const seg_index U; 
      const bool gaps_as_symbols;
      const bool use_pbwt;
      // n[y] is min size ending at column y
      vector<seg_index> n;
      // Traceback
      vector<seg_index> back;
      // Min query
      vector<vector<seg_index>> rings;
      vector<RMaxQTree> rmqtrees;
      vector<RMQueue> rmqueues;
      vector<i_type> keys;
      // Returns the index and size for the min size of M[height][left..right]
      pair<seg_index, seg_index> min_query(const seg_index height, const seg_index y, const seg_index left, const seg_index right) {
        const seg_index l = (max(left - 1, seg_index(0))) % U;
        const seg_index r = (max(right - 1, seg_index(0))) % U;
        seg_index minq = left == 0 ? 0 : numeric_limits<seg_index>::max();
        seg_index x = 0;
        // Naive ring buffer
        if  constexpr (min_range == RING) {
          if (l <= r) {
            for (seg_index j = l; j <= r; j++) {
              if (minq > rings[height - 1][j]) {
                minq = rings[height - 1][j];
                x = left + j - l;
              }
            }
          } else {
            for (seg_index j = l; j < U; j++) {
              if (minq > rings[height - 1][j]) {
                minq = rings[height - 1][j];
                x = left + j - l;
              }
            }
            for (seg_index j = 0; j <= r; j++) {
              if (minq > rings[height - 1][j]) {
                minq = rings[height - 1][j];
                x = left + j + U - l;
              }
            }
          }
        }
        if constexpr (min_range == TREE) {
          if (l <= r) {
            auto [x1, minq1] = rmqtrees[height - 1].query(l + 1, r + 1);
            minq1 *= -1;
            x1--;
            if (minq1 < minq) {
              minq = minq1;
              x = left + x1 - l;
            }
          } else {
            auto [x1, minq1] = rmqtrees[height - 1].query(l + 1, U);
            minq1 *= -1;
            x1--;
            auto [x2, minq2] = rmqtrees[height - 1].query(1, r + 1);
            minq2 *= -1;
            x2--;
            if (minq1 < minq) {
              minq = minq1;
              x = left + x1 - l;
            }
            if (minq2 < minq) {
              minq = minq2;
              x = left + x2 + U - l;
            }
          }
        }
        if constexpr (min_range == QUEUE) {
          long long queue_col = std::max(0LL, y - U);
          long long x_q = rmqueues[height - 1].query(left - queue_col, right - queue_col);
          minq = rmqueues[height - 1].get(x_q);
          x = x_q + queue_col;
        }
        return {x, minq};
      }
       // Updates the min query data structure for row h, column y
      void update(const seg_index h, const seg_index y, const seg_index value) {
        // Naive ring buffer
        if constexpr (min_range == RING) {
          seg_index j = (y - 1) % U;
          rings[h - 1][j] = value;
        }
        if constexpr (min_range == TREE) {
          seg_index j = (y - 1) % U;
          rmqtrees[h - 1].update(j + 1, j + 1, -value);
        }
        if constexpr (min_range == QUEUE) {
          rmqueues[h - 1].push(value);
          if (y >= U)
            rmqueues[h - 1].pop();
        }
      }
      // Calculate the min-L-size for the next column (y is 1-index) given meaningful left extensions L_yy
      void next(const seg_index y, const vector<pair<seg_index, seg_index>> *L_yy) {
        for (size_t j = 0; j + 1 < L_yy->size(); j++) {
          seg_index height = (*L_yy)[j].second;
          seg_index left = (*L_yy)[j + 1].first;
          seg_index right = (*L_yy)[j].first - 1;
          seg_index x, s;

          tie(x, s) = (left == right) 
            ? pair<seg_index, seg_index>{ left, n[left] - height * left } 
            : min_query(height, y, left, right);  
          if (n[y] > height * y + s && height * y + s > 0) {
            n[y] = height * y + s;
            back[y] = x;
          }
        }
        for (seg_index height = 1; height <= r; height++) {
          update(height, y, n[y] - height * y);
        }
      }
      // Calculate the min-L-size for the next column (y is 1-index) using susffix tree algo
      void next_gaps(const seg_index y) {
        set<string> reverse_unique_chunk;
        for (seg_index i = 0; i < r; i++) {
          string s = idx.msa_substr(i, y - min(U, y), min(U, y));
          if (!gaps_as_symbols)
            s.erase(remove(s.begin(), s.end(), '-'), s.end());
          reverse(s.begin(), s.end());
          reverse_unique_chunk.insert(move(s));
        }

        trie::trie T(reverse_unique_chunk);
        reverse_unique_chunk.clear();
        set<string> reverse_chunk;
        for (seg_index i = 0; i < r; i++) {
          string s = idx.msa_substr(i, y - min(U, y), min(U, y));
          reverse(s.begin(), s.end());
          reverse_chunk.insert(s);
        }

        vector<unsigned long long> counts(T.nodes() + 1, 0);
        vector<unsigned long long> node_length(T.nodes() + 1, 0); // length from the root
        counts[T.nodes()] = reverse_chunk.size(); // current count of the root
        vector<unsigned long long> active_node(reverse_chunk.size(), T.nodes());
        for (seg_index len = 1; len <= min(U, y); len++) {
          auto it = reverse_chunk.begin();
          for (seg_index i = 0; i < reverse_chunk.size(); ++i, ++it) {
            if ((*it)[len - 1] != '-') {
              counts[active_node[i]] -= 1;
              auto active_length = node_length[active_node[i]];
              if (active_node[i] == T.nodes()) {
                active_node[i] = T.child(T.root(), (*it)[len - 1]);
              } else {
                active_node[i] = T.child(active_node[i], (*it)[len - 1]);
              }
              assert(active_node[i] != trie::trie::null);
              counts[active_node[i]] += 1;
              node_length[active_node[i]] = active_length + 1;
            }
          }
          assert(it == reverse_chunk.end());    
          // Calculate the total length of the active nodes (segmentation size)
          // Update min-size arrays n and back
          if (len >= L) {
            long long segment_size = 0;
            for (seg_index i = 0; i <= T.nodes(); i++) {
              if (counts[i] > 0) {
                // empty string has size 1
                segment_size += max(1ULL, node_length[i]);
              }
            }
            if (n[y] > n[y - len] + segment_size && n[y - len] + segment_size > 0) {
              n[y] = n[y - len] + segment_size;
              back[y] = y - len;
            }
          }
        }
      }
    public:
      minsize(
        msa_t &idx,
        const seg_index r, 
        const seg_index c, 
        const seg_index L, 
        const seg_index U, 
        const bool gaps_as_symbols,
        const bool use_pbwt
      ): idx(idx), r(r), c(c), L(L), U(U), gaps_as_symbols(gaps_as_symbols), use_pbwt(use_pbwt), 
        n(c + 1, numeric_limits<seg_index>::max()), back(c + 1, 0) {
          n[0] = 0;
          if constexpr (min_range == RING) {
            rings = vector<vector<seg_index>>(r, vector<seg_index>(U, 0));
          } 
          if constexpr (min_range == TREE) {
            rmqtrees.resize(r);
            keys.resize(U + 1);
            for (i_type i = 0; i <= U; ++i)
              keys[i] = i;
            for(auto &rmq: rmqtrees) {
              rmq.fillRMaxQTree(keys.data(), U + 1);
              rmq.update(0, 0, 0);
            }
          }
          if constexpr (min_range == QUEUE) {
            rmqueues = vector<RMQueue>(r, RMQueue(U + 1));
            for(auto &rmq: rmqueues) {
              rmq.push(0);
            }
          }
        }
      /* find the minimum-size segmentation of MSA[1..r][1..c] (indexed by
       *   idx) respecting lower bound L, upper bound U, gaps_as_symbols strategy,
       *   optionally given the meaningful left extensions L_y
       * returns: pair {size, segmentation}, with size =
       *   numeric_limits<seg_index>::max() if no segmentation exists */
      pair<seg_index, vector<pair<seg_index, seg_index>>> segment(
        const vector<vector<pair<seg_index, seg_index>>> &L_y = {}
      ) {
        if (gaps_as_symbols) {
          // algorithm for gaps as symbols from the paper
          for (i_type y = 1; y <= c; y++) {
            // compute L_y if it was not given in input
            const vector<pair<seg_index, seg_index>> *L_yy;
            vector<pair<seg_index, seg_index>> L_yy_on_the_fly;
            if (L_y.size() > 0) {
              L_yy = &(*(L_y.begin() + y));
            } else {
              L_yy_on_the_fly = algo::compute_meaningful_extensions(idx, r, c, L, U, y, gaps_as_symbols, use_pbwt);
              L_yy = &L_yy_on_the_fly;
            }
            // Algorithm for column y
            next(y, L_yy);
          }
        } else {
          // suffix trie algorithm for gaps
          for (i_type y = 1; y <= c; y++) {
            next_gaps(y);
          }
        }

        // trace back
        vector<pair<seg_index, seg_index>> segments;
        for (i_type pos = c; pos > 0; pos = back[pos]) {
          segments.emplace_back(back[pos] + 1, pos);
        }
        reverse(segments.begin(), segments.end());

        return {n[c], segments};
      }
  
  };
}

#endif
