#ifndef BLOCK_GRAPH_HPP
#define BLOCK_GRAPH_HPP
#include <unordered_map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <tuple>

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::ifstream, std::ofstream;
using std::cerr, std::endl;
using std::pair, std::tuple;
using std::max;

// code adapted from https://github.com/algbio/founderblockgraphs/tree/rewrite
namespace eds::block_graph {
    const char GAP_CHARACTER = '-';
    typedef long long seg_index;
    typedef std::vector<long long>::size_type seg_size_t;
    const long long SEG_INDEX_MAX = std::numeric_limits<seg_index>::max();
    typedef vector<pair<seg_index,seg_index>> segmentation;

    /* block graph with nodes partitioned into blocks and arbitrary edges, actually (so a layered DAG) */
    struct block_graph {
        vector<unordered_map<string,unsigned long>> blocks; // each block is a sequence->node id map
        unordered_map<unsigned long,unsigned long> node_to_block; // node id -> its block
        unordered_map<unsigned long,unordered_set<unsigned long>> adjacency_lists; // node id -> out-neighbors
    };

    /* requires: segmentation S is sorted vector of pairs starting at (1,x) and ending at (y,n)
     * returns: elastic block graph (or a layered DAG if a segments contains the empty string)
     * notes: MSA is streamed from disk, graph (no paths) is kept in memory */
    tuple<block_graph,seg_index,seg_index> segment_msa(const string &msa_path, const long long n, const segmentation &S) {
        assert(S.at(0).first == 1 and S.back().second == n);
#ifdef BLOCK_GRAPH_HPP_DEBUG
        cerr << "DEBUG: segmentation segment starts are ";
        for (auto &[s, _] : S) cerr << " " << s;
        cerr << endl;
#endif

        unsigned long nodes = 0;
        vector<unordered_map<string,unsigned long>> blocks(S.size());
        unordered_map<unsigned long,unsigned long> node_to_block;
        unordered_map<unsigned long,unordered_set<unsigned long>> adjacency_lists;
        seg_index card = 0, size = 0; // gap-aware size

        ifstream msa_if(msa_path);
        string line = "", sequence = "";
        while (getline(msa_if, line)) {
            if (line.length() == 0) continue;
            if (line[0] == '>' and sequence.length() > 0) {
                assert(sequence.size() == n);

                seg_index prev = SEG_INDEX_MAX;
                for (seg_size_t i = 0; i < S.size(); i++) {
                    assert(S[i].first <= S[i].second);

                    string label = sequence.substr(S[i].first - 1, S[i].second - S[i].first + 1);
                    label.erase(std::remove(label.begin(), label.end(), GAP_CHARACTER), label.end()); // remove gaps

                    //if (label == "") // we allow gaps
                    //    continue;

                    if (blocks[i].find(label) != blocks[i].end()) {
                        // node label in this block is already present
                        const unsigned long id = blocks[i][label];
                        if (prev != SEG_INDEX_MAX) {
                            assert(id != 0);
                            adjacency_lists[prev].insert(id);
                        }
                        prev = id;
                    } else {
                        // new node
                        const unsigned long newid = nodes++;
                        blocks[i].insert({ label, newid });
                        card += 1;
                        size += max(label.size(), 1LU);
                        node_to_block.insert({ newid, i });
                        adjacency_lists.insert({ newid, unordered_set<unsigned long>() });
                        if (prev != SEG_INDEX_MAX) {
                            assert(newid != 0);
                            adjacency_lists[prev].insert(newid);
                        }
                        prev = newid;
                    }
                }

                sequence = "";
            } else if (line[0] == '>' and sequence.length() == 0) {
                // do nothing
            } else {
                sequence += line;
            }
        }
        if (sequence.size() > 0) {
            seg_index prev = SEG_INDEX_MAX;
            for (seg_size_t i = 0; i < S.size(); i++) {
                assert(S[i].first <= S[i].second);

                string label = sequence.substr(S[i].first - 1, S[i].second - S[i].first + 1);
                label.erase(std::remove(label.begin(), label.end(), GAP_CHARACTER), label.end()); // remove gaps

                //if (label == "") // we allow gaps
                //    continue;

                if (blocks[i].find(label) != blocks[i].end()) {
                    // node label in this block is already present
                    const unsigned long id = blocks[i][label];
                    if (prev != SEG_INDEX_MAX) {
                        assert(id != 0);
                        adjacency_lists[prev].insert(id);
                    }
                    prev = id;
                } else {
                    // new node
                    const unsigned long newid = nodes++;
                    blocks[i].insert({ label, newid });
                    card += 1;
                    size += max(label.size(), 1LU);
                    node_to_block.insert({ newid, i });
                    adjacency_lists.insert({ newid, unordered_set<unsigned long>() });
                    if (prev != SEG_INDEX_MAX) {
                        assert(newid != 0);
                        adjacency_lists[prev].insert(newid);
                    }
                    prev = newid;
                }
            }

            sequence = "";
        }

#ifdef BLOCK_GRAPH_HPP_DEBUG
        cerr << "DEBUG: blocks are ";
        for (auto &b : blocks) {
            cerr << "{";
            for (auto &[label,node] : b) {
                cerr << " " << node << ":" << label;
            }
            cerr << " }";
        }
        cerr << endl;
#endif

        return { block_graph({ std::move(blocks), std::move(node_to_block), std::move(adjacency_lists) }), card, size };
    }

    void output_msa_info(const long long m, const long long n, ofstream &out) {
        out << "M\t" << m << "\t" << n << "\n";
    }
    void output_segmentation(const segmentation &S, ofstream &out) {
        // 0-indexed to 1-indexed, only starting cols (see xGFAspec.md)
        out << "X";
        for (seg_size_t i = 0; i < S.size() - 1; i++)
            out << "\t" << S[i].first;
        out << "\n";
    }
    void output_block_info(const block_graph &g, ofstream &out) {
        out << "B";
        for (auto &b : g.blocks)
            out << "\t" << b.size();
        out << "\n";
    }
    /* TODO: rename vertices? */
    void output_block_graph(const block_graph &g, ofstream &out) {
        for (auto &b : g.blocks) {
            for (auto &[label, node] : b) {
                out << "S\t" << node << "\t" << ((label == "") ? "*" : label) << "\n";
                for (const auto &outneighbor : g.adjacency_lists.at(node)) {
                    out << "L\t" << node << "\t+\t" << outneighbor << "\t+\t0M" << "\n";
                }
            }
        }
    }
    void output_eds(const block_graph &g, ofstream &out) {
        for (auto &b : g.blocks) {
            out << "{";
            bool first = true;
            for (auto &[label, _] : b) {
                out << ((first) ? "" : ",") << label;
                first = false;
            }
            out << "}";
        }
	out << "\n";
    }
} // namespace eds::block_graph
#endif // BLOCK_GRAPH_HPP
