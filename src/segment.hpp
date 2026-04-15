#ifndef BLOCK_GRAPH_HPP
#define BLOCK_GRAPH_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <tuple>
#include <algorithm>

#include "msa_chunker.hpp"

typedef msa_chunker::msa_chunker msa_t;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::ifstream, std::ostream, std::ofstream;
using std::cerr, std::endl;
using std::pair, std::tuple;
using std::max;

namespace segment {
  const char GAP_CHARACTER = '-';
  typedef long long seg_index;
  typedef std::vector<long long>::size_type seg_size_t;
  const long long SEG_INDEX_MAX = std::numeric_limits<seg_index>::max();
  typedef vector<pair<seg_index,seg_index>> segmentation;

  void output_msa_info(const long long m, const long long n, ostream *out) {
    *out << "M\t" << m << "\t" << n << "\n";
  }
  void output_segmentation(const segmentation &S, ostream *out) {
    // 0-indexed to 1-indexed, only starting cols (see xGFAspec.md)
    *out << "X";
    for (seg_size_t i = 0; i < S.size() - 1; i++)
      *out << "\t" << S[i].first;
    *out << "\n";
  }
  void output_block_info(const vector<unsigned long> &block_sizes, ostream *out) {
    *out << "B";
    for (auto size : block_sizes)
      *out << "\t" << size;
    *out << "\n";
  }

  /* requires: segmentation S is sorted vector of pairs starting with [1,x] and ending with [y,c]
   * writes: gfa to out_file
   * returns: cardinality and size of the eds
   * notes: MSA is streamed in chunks from idx, output is streamed */
  pair<seg_index,seg_index> segment_stream_gfa(
      msa_t &idx,
      seg_index rows,
      seg_index cols,
      const vector<pair<seg_index, seg_index>> &S,
      ostream *out
      ) {
    unordered_map<string,unsigned long> block, last_block;
    unordered_map<unsigned long,unordered_set<unsigned long>> adjacency_lists;
    vector<unsigned long> block_sizes;

    seg_index nodes = 0, card = 0, size = 0; // gap-aware size

    output_msa_info(rows, cols, out);
    output_segmentation(S, out);
    block_sizes.reserve(S.size());

    vector<seg_index> prev(rows, SEG_INDEX_MAX);
    for (seg_size_t i = 0; i < S.size(); i++) {
      assert(S[i].first <= S[i].second);
      block.clear();
      adjacency_lists.clear();

      for (seg_size_t r = 0; r < rows; r++) {
        string label = idx.msa_substr(r, S[i].first - 1, S[i].second - S[i].first + 1);
        label.erase(std::remove(label.begin(), label.end(), GAP_CHARACTER), label.end()); // remove gaps

        if (block.find(label) != block.end()) {
          // node label in this block is already present
          const unsigned long id = block[label];
          if (prev[r] != SEG_INDEX_MAX) {
            assert(id != 0);
            adjacency_lists[prev[r]].insert(id);
          }
          prev[r] = id;
        } else {
          // new node
          const unsigned long newid = nodes++;
          block.insert({ label, newid });
          card += 1;
          size += max(label.size(), 1LU);
          adjacency_lists.insert({ newid, unordered_set<unsigned long>() });
          if (prev[r] != SEG_INDEX_MAX) {
            assert(newid != 0);
            adjacency_lists[prev[r]].insert(newid);
          }
          prev[r] = newid;
        }
      }

      // print node labels and edges to this block
      for (auto &[label, node] : block) {
        *out << "S\t" << node << "\t" << ((label == "") ? "*" : label) << "\n";
      }
      for (auto &[inneighbor, list] : adjacency_lists) {
        for (auto &outneighbor : list) {
          *out << "L\t" << inneighbor << "\t+\t" << outneighbor << "\t+\t0M" << "\n";
        }
      }
      block_sizes.push_back(block.size());
      last_block = std::move(block);
    }

    output_block_info(block_sizes, out);
    return { card, size };
  }

  pair<seg_index,seg_index> segment_stream_eds(
      msa_t &idx,
      seg_index rows,
      seg_index cols,
      const vector<pair<seg_index, seg_index>> &S,
      ostream *out
      ) {
    seg_index nodes = 0, card = 0, size = 0; // gap-aware size
    std::set<string> block;

    for (seg_size_t i = 0; i < S.size(); i++) {
      assert(S[i].first <= S[i].second);
      block.clear();

      for (seg_size_t r = 0; r < rows; r++) {
        string label = idx.msa_substr(r, S[i].first - 1, S[i].second - S[i].first + 1);
        label.erase(std::remove(label.begin(), label.end(), GAP_CHARACTER), label.end()); // remove gaps
        block.insert(label);
      }

      // print block
      card += block.size();
      *out << "{";
      bool first = true;
      for (auto &label : block) {
        size += max(label.size(), 1LU);
        *out << ((first) ? "" : ",") << label;
        first = false;
      }
      *out << "}";
    }
    *out << "\n";
    return { card, size };
  }

  pair<seg_index,seg_index> segment_stream_no_output(
      msa_t &idx,
      seg_index rows,
      seg_index cols,
      const vector<pair<seg_index, seg_index>> &S
      ) {
    seg_index nodes = 0, card = 0, size = 0; // gap-aware size
    std::set<string> block;

    for (seg_size_t i = 0; i < S.size(); i++) {
      assert(S[i].first <= S[i].second);
      block.clear();

      for (seg_size_t r = 0; r < rows; r++) {
        string label = idx.msa_substr(r, S[i].first - 1, S[i].second - S[i].first + 1);
        label.erase(std::remove(label.begin(), label.end(), GAP_CHARACTER), label.end()); // remove gaps
        block.insert(label);
      }

      card += block.size();
      for (auto &label : block) {
        size += max(label.size(), 1LU);
      }
    }
    return { card, size };
  }
} // namespace segment
#endif // BLOCK_GRAPH_HPP
