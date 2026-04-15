#ifndef TRIE_HPP
#define TRIE_HPP

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <iterator>
#include <set>
#include <tuple>
#include <limits>

#include "sdsl/bp_support.hpp"
#include "sdsl/rank_support.hpp"

using std::vector;
using std::string;
using std::set;
using std::tuple;
using std::cerr, std::endl;
using std::reference_wrapper;

//#define TRIE_DEBUG

/*
 * static string trie allowing construction, navigation, and node comparisons
 */
namespace trie {
  // non-root node indexes are in the range [0,nodes())
  class trie {
    private:
      typedef sdsl::bp_support_sada<> bp_support_t;
      const unsigned long long root_i = std::numeric_limits<unsigned long long>::max();
      unsigned long long _nodes, marked_nodes;
      sdsl::bit_vector topology;
      bp_support_t bp_topology;
      string label_to_parent;

      inline unsigned long long lcp(const string &s1, const string &s2) {
        unsigned long long res = 0;
        while (res < s1.length() and res < s2.length() and s1[res] == s2[res])
          res += 1;
        return res;
      }

    public:
      static const unsigned long long null = std::numeric_limits<unsigned long long>::max();
      trie() = delete;

      trie(const set<string> &in) {
        _nodes = 0;
        marked_nodes = 0;

        if (in.size() == 0) {
          return;
        }

        // find LCP values
        vector<unsigned long long> lcp_lengths(in.size());
        lcp_lengths[0] = 0;
        _nodes += in.begin()->length();
        for (auto s = ++(in.begin()), prev = in.begin(); s != in.end(); ++s, ++prev) {
          const unsigned long long l = lcp(*prev, *s);
#ifdef TRIE_DEBUG
          cerr << "DEBUG: string " << *s << " LCP: " << l << std::endl;
#endif
          lcp_lengths[++marked_nodes] = l;
          _nodes += s->length() - l;
        }
#ifdef TRIE_DEBUG
        cerr << "DEBUG: lcp array is";
        for (auto i : lcp_lengths)
          cerr << " " << i;
        cerr << endl;
#endif

        // find balanced parenthesis representation
        topology = sdsl::bit_vector(2 * _nodes, 0);
        label_to_parent.resize(_nodes);
        auto s = in.begin();
        std::size_t n = 0, nn = 0, length = s->length();
        while (n < length) {
          topology[n] = 1; // open the nodes of the first string
          label_to_parent[nn] = (*s)[n];
          n += 1;
          nn += 1;
        }
        std::size_t prevlength = length;
        ++s;
        for (std::size_t i = 1; i < lcp_lengths.size(); ++i, ++s) {
          if (prevlength != lcp_lengths[i]) { // string is not prefix of previous one
            for (std::size_t j = 0; j < prevlength - lcp_lengths[i]; ++j) {
              topology[n++] = 0; // close the exhausted nodes of the previous string
            }
          }

          length = s->length();
          const std::size_t i_lcp_length = lcp_lengths[i];
          for (std::size_t j = 0; j < length - i_lcp_length; ++j) {
            topology[n++] = 1; // open the new nodes of this string
            label_to_parent[nn++] = (*s)[i_lcp_length + j];
          }

          prevlength = length;
        }
        // currently-opened nodes are already 0s in topology
        assert(s == in.end());
#ifdef TRIE_DEBUG
        cerr << "DEBUG: topology array is " << topology << endl;
        cerr << "DEBUG: label_to_parent is " << label_to_parent << endl;
#endif
        bp_topology = bp_support_t(&topology);
      }

      unsigned long long nodes() const {
        return _nodes;
      }

      unsigned long long root() const {
        return root_i;
      }

      unsigned long long child(unsigned long long n, char c) {
        const unsigned long long n_i = ((n == root_i) ? root_i : bp_topology.select(n + 1)); // open parenthesis index
        const unsigned long long upper_bound_i = ((n == root_i) ? topology.size() : bp_topology.find_close(n_i)) ; // close parenthesis index
        const unsigned long long first_child_i = ((n == root_i) ? 0 : n_i + 1);
        for (unsigned long long i = first_child_i; i < upper_bound_i; i = bp_topology.find_close(i) + 1) {
          const unsigned long long nn = bp_topology.rank(i) - 1;
          if (label_to_parent[nn] == c)
            return nn;
        }
        return null; // no child was labeled c
      }
  };
} // namespace trie
#endif // TRIE_HPP
