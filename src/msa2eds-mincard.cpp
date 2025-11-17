#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <algorithm>
#include <climits>
#include <unordered_set>
#include <chrono>
#include <limits>

#include "RMaxQTree.h"
#include "block_graph.hpp"

using namespace std::chrono;
using namespace std;
using eds::block_graph::block_graph, eds::block_graph::segment_msa, eds::block_graph::output_msa_info, eds::block_graph::output_segmentation, eds::block_graph::output_block_info, eds::block_graph::output_block_graph;

bool verbose = false;
typedef eds::block_graph::seg_index seg_index;
typedef long long int key_type;

// Reads sequences from a FASTA file
vector<string> read_fasta(const string& filename) {
    ifstream in(filename);
    vector<string> sequences;
    string line, current;

    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current.empty()) {
                sequences.push_back(current);
                current.clear();
            }
        } else {
            current += line;
        }
    }
    if (!current.empty())
        sequences.push_back(current);

    return sequences;
}

vector<vector<pair<seg_index, seg_index>>> compute_meaningful_extensions(
    const vector<string>& msa, seg_index L, seg_index U)
{
    seg_index r = msa.size();
    seg_index c = msa[0].size();

    vector<vector<pair<seg_index, seg_index>>> L_y(c + 1);  // 1-based indexing

    for (seg_index y = 1; y <= c; ++y) {
        if (y < L) {
            continue;  // No extension possible
        }

        vector<pair<seg_index, seg_index>> current;

        // Enforce ℓ_{y,1} = y - L + 1 down to ℓ_{y,d_y} > y - U
        seg_index prev_height = -1;
        for (seg_index len = L; len <= U && y - len + 1 >= 1; ++len) {
            seg_index start = y - len + 1;
            unordered_set<string> unique_strings;

            for (seg_index i = 0; i < r; ++i) {
                string s = msa[i].substr(start - 1, len);
                s.erase(remove(s.begin(), s.end(), '-'), s.end());  // remove gaps
                unique_strings.insert(s);
            }

            seg_index height = unique_strings.size();
            if (height != prev_height) {
                current.emplace_back(start, height);
                prev_height = height;
            }
        }

        // Add dummy ℓ_{y,d_y+1} = max(0, y - U)
        current.emplace_back(max((seg_index)0, y - U), -1);
        L_y[y] = current;
    }

    return L_y;
}

vector<bool> compute_perfect_columns(
    const vector<string>& msa) {
    seg_index r = msa.size();
    seg_index c = msa[0].size();
    assert(r > 0);

    vector<bool> perfect_columns(c + 1, true); // 1-indexed
    for (seg_index y = 1; y <= c; ++y) {
        const char consensus = msa[0][y-1];
        for (seg_index i = 2; i <= r; ++i) {
            if (msa[i-1][y-1] != consensus) {
                perfect_columns[y] = false;
                break;
            }
        }
    }
    return perfect_columns;
}

const vector<bool> perfect_columns_dummy = {};
pair<seg_index, vector<pair<seg_index, seg_index>>> segment_with_rmq(
    const vector<vector<pair<seg_index, seg_index>>>& L_y, seg_index c, const vector<bool> &perfect_columns = perfect_columns_dummy)
{
    const bool allow_perfect_segments = (perfect_columns.size() > 0);
    vector<seg_index> m(c + 1, numeric_limits<seg_index>::max());      // m[y] is the DP value: minimal number of strings
    vector<seg_index> mneg(c + 1, numeric_limits<seg_index>::min());  // store -m[y] for max-query simulation
    vector<seg_index> back(c + 1, -1);    // traceback
    seg_index perfect_back = -1, perfect_m = numeric_limits<seg_index>::max();
    if (allow_perfect_segments and perfect_columns[0]) {
        perfect_m = 0;
        perfect_back = 0;
    }

    m[0] = 0;
    mneg[0] = 0;

    // Initial fill of keys = 0..c
    vector<key_type> keys(c + 1);
    for (key_type i = 0; i <= c; ++i) keys[i] = i;

    // Initialize RMaxQTree with negated m-values
    RMaxQTree rmq;
    rmq.fillRMaxQTree(keys.data(), c + 1);
    rmq.update(0, 0, 0);  // set index 0 with mneg[0] = 0

    for (key_type y = 1; y <= c; ++y) {
        m[y] = numeric_limits<key_type>::max();

        const auto& L = L_y[y];

        // optimal solution using L_y
        for (size_t j = 0; j + 1 < L.size(); ++j) {
            key_type l = L[j + 1].first;
            key_type r = L[j].first - 1;
            if (l > r) continue;

            // query returns pair (index, value), but value is -m[index]
            auto [x, neg_mx] = rmq.query(l, r);
            key_type candidate = L[j].second + m[x];

            if (candidate < m[y]) {
                m[y] = candidate;
                back[y] = x;
            }
        }

        if (allow_perfect_segments and perfect_columns[y]) {
            if (perfect_m < numeric_limits<seg_index>::max() and perfect_m + 1 <= m[y]) {
                m[y] = perfect_m + 1;
                back[y] = perfect_back;
            }
        }

        mneg[y] = -m[y];
        rmq.update(y, y, mneg[y]);

        // optional
        if (allow_perfect_segments and y < c and perfect_columns[y+1]) {
            if (m[y] < perfect_m) {
                perfect_m = m[y];
                perfect_back = y;
            }
        } else if (allow_perfect_segments and y < c and !perfect_columns[y+1]) {
            perfect_m = numeric_limits<seg_index>::max();
            perfect_back = -1;
        }
    }

    // Traceback
    vector<pair<seg_index, seg_index>> segments;
    for (key_type pos = c; pos > 0; pos = back[pos]) {
        segments.emplace_back(back[pos] + 1, pos);
    }
    reverse(segments.begin(), segments.end());

    return {m[c], segments};
}

// Prseg_index EDS from segmentation
void prseg_index_eds(const vector<string>& msa, const vector<pair<seg_index, seg_index>>& segments, string out_filename = "") {
    std::ofstream outFile;
    if (out_filename.size()==0) 
        cout << "Elastic Degenerate String (EDS):\n";
    else
       outFile.open(out_filename);
    for (const auto& [l, r] : segments) {
        set<string> unique_subs;
        for (const auto& seq : msa) {
            string sub = seq.substr(l - 1, r - l + 1);
            sub.erase(remove(sub.begin(), sub.end(), '-'), sub.end()); // remove gaps
            unique_subs.insert(sub);
        }
        if (out_filename.size()==0) { 
            cout << "{ ";
            for (const auto& s : unique_subs)
                cout << s << " ";
            cout << "}";
        }
        if (outFile) {
            outFile << "{ ";
            for (const auto& s : unique_subs)
                outFile << s << " ";
            outFile << "}";
        }         
    }
    if (out_filename.size()==0) 
        cout << "\n";
    if (outFile)
        outFile.close();
}

// Count the total cardinality of sets
seg_index card_eds(const vector<string>& msa, const vector<pair<seg_index, seg_index>>& segments) {
    seg_index card = 0;
    for (const auto& [l, r] : segments) {
        set<string> unique_subs;
        for (const auto& seq : msa) {
            string sub = seq.substr(l - 1, r - l + 1);
            sub.erase(remove(sub.begin(), sub.end(), '-'), sub.end()); // remove gaps
            unique_subs.insert(sub);
        }
        card += unique_subs.size();
    }
    return card;
}

// Main function
int main(int argc, char* argv[]) {
    string filename = "example.fasta";
    seg_index L = 1;
    seg_index U = 10;
    bool allow_perfect_segments = false;

    if (argc<=1) {
      cout << "Syntax: " << string(argv[0]) << " msa.fasta segment-length-upper-bound (default " << U << ") allow-perfect-segments (default 0) verbose (default 0)" << endl;
      return 0;
    }

    filename = string(argv[1]);
    if (argc>2)
      U = atoi(argv[2]);
    if (argc>3)
      allow_perfect_segments = atoi(argv[3]) > 0;
    if (argc>4)
      verbose = atoi(argv[4]);
    auto msa = read_fasta(filename);
    if (msa.empty()) {
        cerr << "MSA file is empty or not found.\n";
        return 1;
    }
    auto start_pre = high_resolution_clock::now();
    auto L_y = compute_meaningful_extensions(msa, L, U);
    auto stop_pre = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop_pre-start_pre);
    cout << "Preprocessing took " << duration.count() << " milliseconds" << endl;
    if (verbose) {
        cout << "Meaningful left extensions and heights:\n";
        for (seg_index y = 1; y < L_y.size(); ++y) {
            if (L_y[y].empty()) continue;
            cout << "y = " << y << ":\n";
            for (size_t j = 0; j < L_y[y].size(); ++j) {
                cout << "  ℓ: " << L_y[y][j].first << "   h: " << L_y[y][j].second << "\n";
            }
        }
    }

    vector<bool> perfect_columns = {};
    if (allow_perfect_segments) {
            perfect_columns = compute_perfect_columns(msa);
    }
    auto start_dp = high_resolution_clock::now();
    auto [cost, segments] = segment_with_rmq(L_y, msa[0].size(), perfect_columns);
    auto stop_dp = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_dp-start_dp);
    cout << "DP took " << duration.count() << " milliseconds" << endl;

    cout << "Minimum segmentation cardinality: " << cost << "\n";
    if (verbose) {
       cout << "Segments:\n";
       for (auto [l, r] : segments)
           cout << "[" << l << "," << r << "] ";
       cout << "\n";

       prseg_index_eds(msa, segments);
    }
    //prseg_index_eds(msa,segments,filename + ".eds.txt");
    block_graph eds = segment_msa(filename, msa[0].size(), segments);
    ofstream out(filename + ".gfa");
    output_msa_info(msa.size(), msa[0].size(), out);
    output_segmentation(segments, out);
    output_block_info(eds, out);
    output_block_graph(eds, out);
    cout << "Minimum segmentation cardinality after gap removal: " << card_eds(msa, segments) << "\n";
    return 0;
}
