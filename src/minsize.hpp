#ifndef MINSIZE_HPP
#define MINSIZE_HPP

#include "segment.hpp"
#include "msa_chunker.hpp"
#include "algo.hpp"

using std::vector, std::tie;

namespace minsize {
  typedef msa_chunker::msa_chunker msa_t;
  typedef segment::seg_index seg_index;

  // Class for solving the min-L-size segmentation problem 
  class minsize {
      msa_t &idx;
      const seg_index r;
      const seg_index c; 
      const seg_index L; 
      const seg_index U; 
      const bool gaps_as_symbols;
      const bool use_pbwt;
      // Backtracking
      vector<seg_index> back;
      seg_index n_y;
      // Min query
      vector<vector<seg_index>> m;
      // Returns the index and size for the min size of M[height][left..right]
      pair<seg_index, seg_index> min_query(const seg_index height, const seg_index left, const seg_index right) {
        const seg_index l = (max(left - 1, seg_index(0))) % U;
        const seg_index r = (right - 1) % U;
        seg_index minq = left == 0 ? 0 : numeric_limits<seg_index>::max();
        seg_index x = 0;
        // Naive ring buffer
        if (l <= r) {
          for (seg_index j = l; j <= r; j++) {
            if (minq > m[height - 1][j]) {
              minq = m[height - 1][j];
              x = left + j - l;
            }
          }
        } else {
          for (seg_index j = l; j < U; j++) {
            if (minq > m[height - 1][j]) {
              minq = m[height - 1][j];
              x = left + j - l;
            }
           }
          for (seg_index j = 0; j <= r; j++) {
            if (minq > m[height - 1][j]) {
              minq = m[height - 1][j];
              x = left + j + U - l - 1;
            }
          }
        }
        return {x, minq};
      }
       // Updates the min query data structure for row h, column y
      void update(const seg_index h, const seg_index y, const seg_index value) {
        // Naive ring buffer
        seg_index j = (y - 1) % U;
        m[h - 1][j] = value;
      }
      // Calculate the min-L-size for the next column (y is 1-index) given meaningful left extensions L_yy
      void next(const seg_index y, const vector<pair<seg_index, seg_index>> *L_yy) {
        n_y = numeric_limits<seg_index>::max();
        for (size_t j = 0; j + 1 < L_yy->size(); j++) {
          seg_index height = (*L_yy)[j].second;
          seg_index left = (*L_yy)[j + 1].first;
          seg_index right = (*L_yy)[j].first - 1;
          seg_index x, s;
          tie(x, s) = min_query(height, left, right); 
          // if (y < 15)
          //   cerr << "y = " << y << " h = " << height << " l = " << left << " r = " << right << " x = " << x << " s = " << s << " n_y = " << n_y << " small = " <<  height * y + s << endl;
          if (n_y > height * y + s && height * y + s > 0) {
            n_y = height * y + s;
            back[y] = x;
          }
        }
        for (seg_index height = 1; height <= r; height++) {
          update(height, y, n_y - height * y);
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
        n_y(numeric_limits<seg_index>::max()), back(c + 1, 0),
        m(r, vector<seg_index>(U, 0)) {}
      /* find the minimum-size segmentation of MSA[1..r][1..c] (indexed by
       *   idx) respecting lower bound L, upper bound U, gaps_as_symbols strategy,
       *   optionally given the meaningful left extensions L_y
       * returns: pair {size, segmentation}, with size =
       *   numeric_limits<seg_index>::max() if no segmentation exists */
      pair<seg_index, vector<pair<seg_index, seg_index>>> segment(
        const vector<vector<pair<seg_index, seg_index>>> &L_y = {}
      ) {
        // Algorithm for gaps as symbols
        // TODO: Add trie algo for gaps
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

        // trace back
        vector<pair<seg_index, seg_index>> segments;
        for (i_type pos = c; pos > 0; pos = back[pos]) {
          segments.emplace_back(back[pos] + 1, pos);
        }
        reverse(segments.begin(), segments.end());

        return {n_y, segments};
      }
  
  };
}

#endif
