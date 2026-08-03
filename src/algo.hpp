#ifndef ALGO_HPP
#define ALGO_HPP

#include <vector>
#include <tuple>
#include <set>
#include <climits>

#include "segment.hpp"
#include "trie.hpp"
#include "rmqueue.h"
#include "msa_chunker.hpp"
#include "pbwt.h"

using std::vector;
using std::pair;
using std::set;
using std::numeric_limits;

//#define ALGO_DEBUG

namespace algo {
  typedef msa_chunker::msa_chunker msa_t;
  typedef segment::seg_index seg_index;

  /* compute the meaningful extensions L_y of MSA[1..r][1..c] indexed by idx */
  vector<pair<seg_index, seg_index>> compute_meaningful_extensions(
      msa_t &idx,
      const seg_index r,
      const seg_index c,
      const seg_index L,
      const seg_index U,
      const seg_index y,
      const bool gaps_as_symbols,
      const bool use_pbwt
      ) {
    static pbwt pb = pbwt(r);
    if (use_pbwt) {
      pb.next(idx.get_column(y - 1));
      return pb.meaningful_extensions(L, U);
    }

    vector<pair<seg_index, seg_index>> L_y;  // 1-based indexing
    if (y < L) {
      return L_y; // No extension possible
    }
    
    set<string> reverse_unique_chunk;
    for (seg_index i = 0; i < r; i++) {
      string s = idx.msa_substr(i, y - min(U, y), min(U, y));
      if (!gaps_as_symbols)
        s.erase(remove(s.begin(), s.end(), '-'), s.end());
      reverse(s.begin(), s.end());
      reverse_unique_chunk.insert(move(s));
    }

    if (reverse_unique_chunk.size() > 1) {
      trie::trie T(reverse_unique_chunk);
      reverse_unique_chunk.clear();

      set<string> reverse_chunk;
      for (seg_index i = 0; i < r; i++) {
        string s = idx.msa_substr(i, y - min(U, y), min(U, y));
        reverse(s.begin(), s.end());
        reverse_chunk.insert(s);
      }

      vector<unsigned long long> counts(T.nodes() + 1, 0);
      counts[T.nodes()] = reverse_chunk.size(); // current count of the root
      vector<unsigned long long> active_node(reverse_chunk.size(), T.nodes());
      unsigned long long height = 1; // distinct nodes

      for (seg_index len = 1; len <= min(U, y); ++len) {
        unsigned long long prev_height = height;
        auto it = reverse_chunk.begin();
        for (seg_index i = 0; i < reverse_chunk.size(); ++i, ++it) {
          if (gaps_as_symbols or (*it)[len - 1] != '-') {
            counts[active_node[i]] -= 1;
            if (counts[active_node[i]] == 0) {
              height -= 1;
            }

            if (active_node[i] == T.nodes()) {
              active_node[i] = T.child(T.root(), (*it)[len - 1]);
            } else {
              active_node[i] = T.child(active_node[i], (*it)[len - 1]);
            }
            if (counts[active_node[i]] == 0) {
              height += 1;
            }
            counts[active_node[i]] += 1;
            assert(active_node[i] != trie::trie::null);
          }
        }
        assert(it == reverse_chunk.end());
        if (len == L or (len > L and prev_height != height)) {
          L_y.emplace_back(y - len + 1, height);
        }
      }
    } else {
      if (gaps_as_symbols) {
        L_y.emplace_back(y - L + 1, 1);
      } else {
        set<string> reverse_chunk;
        for (seg_index i = 0; i < r; i++) {
          string s = idx.msa_substr(i, y - min(U, y), min(U, y));
          reverse(s.begin(), s.end());
          reverse_chunk.insert(s);
        }
        if (reverse_chunk.size() == 1) {
          L_y.emplace_back(y - L + 1, 1);
        } else {
          vector<unsigned long long> counts(min(U, y) + 1, 0);
          counts[0] = reverse_chunk.size(); // current count of ε
          vector<unsigned long long> active_node(reverse_chunk.size(), 0);
          unsigned long long height = 1; // distinct nodes

          for (seg_index len = 1; len <= min(U, y); ++len) {
            unsigned long long prev_height = height;
            auto it = reverse_chunk.begin();
            for (seg_index i = 0; i < reverse_chunk.size(); ++i, ++it) {
              if ((*it)[len - 1] != '-') {
                counts[active_node[i]] -= 1;
                if (counts[active_node[i]] == 0) {
                  height -= 1;
                }

                active_node[i] += 1;
                if (counts[active_node[i]] == 0) {
                  height += 1;
                }
                counts[active_node[i]] += 1;
              }
            }
            assert(it == reverse_chunk.end());
            if (len == L or (len > L and prev_height != height)) {
              L_y.emplace_back(y - len + 1, height);
            }
          }
        }
      }
    }

    // Add dummy ℓ_{y,d_y+1} = max(0, y - U)
    L_y.emplace_back(max((seg_index)0, y - U), -1);

    return L_y;
  }

  /* naive(r) version of compute_meaningful_extensions for debugging */
  vector<pair<seg_index, seg_index>> compute_meaningful_extensions_naive(
      msa_t &idx,
      const seg_index r,
      const seg_index c,
      const seg_index L,
      const seg_index U,
      const seg_index y,
      const bool gaps_as_symbols
      ) {
    vector<pair<seg_index, seg_index>> L_y;  // 1-based indexing
    if (y < L) {
      return L_y; // No extension possible
    }

    // Compute ℓ_{y,1}, ℓ_{y,2}, ..., ℓ_{y,d_y} in reverse order
    seg_index prev_height = -1;
    for (seg_index len = min(U, y); len >= L; --len) {
      seg_index start = y - len + 1;
      unordered_set<string> unique_strings;

      for (seg_index i = 0; i < r; ++i) {
        string s = idx.msa_substr(i, start - 1, len);
        if (!gaps_as_symbols)
          s.erase(remove(s.begin(), s.end(), '-'), s.end());
        unique_strings.insert(s);
      }

      seg_index height = unique_strings.size();
      if (prev_height != -1 and height != prev_height) {
        L_y.emplace_back(start - 1, prev_height);
      }
      prev_height = height;
    }
    L_y.emplace_back(y - L + 1, prev_height);
    reverse(L_y.begin(), L_y.end());

    // Add dummy ℓ_{y,d_y+1} = max(0, y - U)
    L_y.emplace_back(max((seg_index)0, y - U), -1);

    return L_y;
  }

  vector<vector<pair<seg_index, seg_index>>> compute_all_meaningful_extensions(
      msa_t &idx,
      const seg_index r,
      const seg_index c,
      const seg_index Lbound,
      const seg_index Ubound,
      const bool gaps_as_symbols,
      const bool use_pbwt
      ) {
    vector<vector<pair<seg_index, seg_index>>> L(c + 1);  // 1-based indexing

    for (seg_index y = 1; y <= c; ++y) {
      L[y] = compute_meaningful_extensions(idx, r, c, Lbound, Ubound, y, gaps_as_symbols, use_pbwt);
    }

    return L;
  }

  /* compute the boolean vector marking perfect columns (i.e. no variation) */
  pair<seg_index,vector<bool>> compute_perfect_columns(
      msa_t &idx, const seg_index r, const seg_index c) {
    seg_index np = 0;
    assert(r > 0);

    vector<bool> perfect_columns(c + 1, true); // 1-indexed
    for (seg_index y = 1; y <= c; ++y) {
      const char consensus = idx.msa_substr(0, y-1, 1)[0];
      for (seg_index i = 1; i < r; ++i) {
        const char i_char = idx.msa_substr(i, y-1, 1)[0];
        if (i_char != consensus) {
          perfect_columns[y] = false;
          np += 1;
          break;
        }
      }
    }
    return {c - np, std::move(perfect_columns)};
  }

  const vector<bool> perfect_columns_dummy = {};
  /* find the minimum-cardinality segmentation of MSA[1..r][1..c] (indexed by
   *   idx) respecting lower bound L, upper bound U, gaps_as_symbols strategy,
   *   optionally given the meaningful left extensions L_y and perfect_column
   *   info (Algorithm 1 in the paper)
   * returns: pair {cardinality, segmentation}, with cardinality =
   *   numeric_limits<seg_index>::max() if no segmentation exists */
  pair<seg_index, vector<pair<seg_index, seg_index>>> segment_with_rmq(
      msa_t &idx,
      const seg_index r,
      const seg_index c,
      const int L,
      const int U,
      const bool gaps_as_symbols,
      const bool use_pbwt,
      const vector<vector<pair<seg_index, seg_index>>> &L_y,
      const vector<bool> &perfect_columns = perfect_columns_dummy
      ) {
    vector<seg_index> m(c + 1); // m[y] is min card ending at column y (1-indexed)
    vector<seg_index> back(c + 1, -1); // traceback for ending segment cols (1-indexed)

    seg_index perfect_first = -1, perfect_m = numeric_limits<seg_index>::max();
    const bool allow_perfect_segments = (perfect_columns.size() > 0);

    RMQueue rmq(U + 1);
    m[0] = 0;
    rmq.push(0);

    if (allow_perfect_segments) {
      perfect_m = numeric_limits<seg_index>::max();
      perfect_first = -1;
    }

    for (seg_index y = 1; y <= c; ++y) {
      m[y] = numeric_limits<seg_index>::max();

      // compute L_y if it was not given in input
      const vector<pair<seg_index, seg_index>> *L_yy;
      vector<pair<seg_index, seg_index>> L_yy_on_the_fly;
      if (L_y.size() > 0) {
        L_yy = &(*(L_y.begin() + y));
      } else {
        L_yy_on_the_fly = compute_meaningful_extensions(idx, r, c, L, U, y, gaps_as_symbols, use_pbwt);
        L_yy = &L_yy_on_the_fly;
      }
#ifdef ALGO_DEBUG
      assert(*L_yy == compute_meaningful_extensions_naive(idx, r, c, L ,U, y, gaps_as_symbols));
#endif

      // optimal solution using L_y
      for (size_t j = 0; j + 1 < L_yy->size(); ++j) {
        seg_index l = (*L_yy)[j + 1].first;
        seg_index r = (*L_yy)[j].first - 1;
        if (l > r) continue;

        seg_index x = l;
        if (L > 1) {
          // m_y is not monotone non-decreasing
          seg_index queue_col = std::max(0LL, y - U);
          x = rmq.query(l - queue_col, r - queue_col);
          auto mx = rmq.get(x);
          
          if (x == -1) continue;
          if (mx == numeric_limits<seg_index>::max()) continue;
          x += queue_col;
        }
        seg_index candidate = (*L_yy)[j].second + m[x];
        assert(candidate >= 0);

        if (candidate < m[y]) {
          m[y] = candidate;
          back[y] = x;
        }
      }

      // optionally use the best in perfect-segment run
      if (allow_perfect_segments and perfect_columns[y]) {
        if (perfect_m < numeric_limits<seg_index>::max() and perfect_m + 1 <= m[y]) {
          m[y] = perfect_m + 1;
          back[y] = perfect_first - 1;
        }
      }

      rmq.push(m[y]);
      if (y >= U)
        rmq.pop();

      // update perfect-segment run
      if (allow_perfect_segments) {
        if (perfect_columns[y] and m[y - 1] < perfect_m) {
          perfect_m = m[y - 1]; // range min query is shifted by 1!
          perfect_first = y;
        } else if (!perfect_columns[y]) {
          perfect_m = numeric_limits<seg_index>::max();
          perfect_first = -1;
        }
      }
    }

    // trace back
    vector<pair<seg_index, seg_index>> segments;
    for (seg_index pos = c; pos > 0; pos = back[pos]) {
      segments.emplace_back(back[pos] + 1, pos);
    }
    reverse(segments.begin(), segments.end());

    return {m[c], segments};
  }
} // namespace algo
#endif // ALGO_HPP
