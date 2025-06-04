#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <algorithm>
#include <climits>
#include "RMaxQTree.h"
#include <unordered_set>
#include <chrono>
using namespace std::chrono;
using namespace std;

const int INF = INT_MAX;
bool verbose = false;

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

vector<vector<pair<int, int>>> compute_meaningful_extensions(
    const vector<string>& msa, int L, int U)
{
    int r = msa.size();
    int c = msa[0].size();

    vector<vector<pair<int, int>>> L_y(c + 1);  // 1-based indexing

    for (int y = 1; y <= c; ++y) {
        if (y < L) {
            continue;  // No extension possible
        }

        vector<pair<int, int>> current;

        // Enforce ℓ_{y,1} = y - L + 1 down to ℓ_{y,d_y} > y - U
        int prev_height = -1;
        for (int len = L; len <= U && y - len + 1 >= 1; ++len) {
            int start = y - len + 1;
            unordered_set<string> unique_strings;

            for (int i = 0; i < r; ++i) {
                string s = msa[i].substr(start - 1, len);
                s.erase(remove(s.begin(), s.end(), '-'), s.end());  // remove gaps
                unique_strings.insert(s);
            }

            int height = unique_strings.size();
            if (height != prev_height) {
                current.emplace_back(start, height);
                prev_height = height;
            }
        }

        // Add dummy ℓ_{y,d_y+1} = max(0, y - U)
        current.emplace_back(max(0, y - U), -1);
        L_y[y] = current;
    }

    return L_y;
}

pair<int, vector<pair<int, int>>> segment_with_rmq(
    const vector<vector<pair<int, int>>>& L_y, int c)
{
    const int INF = 1 << 30;

    vector<int> m(c + 1, INF);      // m[y] is the DP value: minimal number of strings
    vector<int> mneg(c + 1, -INF);  // store -m[y] for max-query simulation
    vector<int> back(c + 1, -1);    // traceback

    m[0] = 0;
    mneg[0] = 0;

    // Initial fill of keys = 0..c
    vector<int> keys(c + 1);
    for (int i = 0; i <= c; ++i) keys[i] = i;

    // Initialize RMaxQTree with negated m-values
    RMaxQTree rmq;
    rmq.fillRMaxQTree(keys.data(), c + 1);
    rmq.update(0, 0, 0);  // set index 0 with mneg[0] = 0

    for (int y = 1; y <= c; ++y) {
        m[y] = INF;

        const auto& L = L_y[y];

        for (size_t j = 0; j + 1 < L.size(); ++j) {
            int l = L[j + 1].first;
            int r = L[j].first - 1;
            if (l > r) continue;

            // query returns pair (index, value), but value is -m[index]
            auto [x, neg_mx] = rmq.query(l, r);
            int candidate = L[j].second + m[x];

            if (candidate < m[y]) {
                m[y] = candidate;
                back[y] = x;
            }
        }

        mneg[y] = -m[y];
        rmq.update(y, y, mneg[y]);
    }

    // Traceback
    vector<pair<int, int>> segments;
    for (int pos = c; pos > 0; pos = back[pos]) {
        segments.emplace_back(back[pos] + 1, pos);
    }
    reverse(segments.begin(), segments.end());

    return {m[c], segments};
}

// Print EDS from segmentation
void print_eds(const vector<string>& msa, const vector<pair<int, int>>& segments, string out_filename = "") {
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
int card_eds(const vector<string>& msa, const vector<pair<int, int>>& segments) {
    int card = 0;
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
    int L = 1;
    int U = 10;

    if (argc<=1) {
      cout << "Syntax: " << string(argv[0]) << " msa.fasta segment-length-upper-bound (default " << U << ") verbose (default 0)" << endl;
      return 0;
    }

    filename = string(argv[1]);
    if (argc>2)
      U = atoi(argv[2]);
    if (argc>3)
      verbose = atoi(argv[3]);
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
        for (int y = 1; y < L_y.size(); ++y) {
            if (L_y[y].empty()) continue;
            cout << "y = " << y << ":\n";
            for (size_t j = 0; j < L_y[y].size(); ++j) {
                cout << "  ℓ: " << L_y[y][j].first << "   h: " << L_y[y][j].second << "\n";
            }
        }
    }
    auto start_dp = high_resolution_clock::now();
    auto [cost, segments] = segment_with_rmq(L_y, msa[0].size());
    auto stop_dp = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop_dp-start_dp);
    cout << "DP took " << duration.count() << " milliseconds" << endl;

    cout << "Minimum segmentation cardinality: " << cost << "\n";
    if (verbose) {
       cout << "Segments:\n";
       for (auto [l, r] : segments)
           cout << "[" << l << "," << r << "] ";
       cout << "\n";

       print_eds(msa, segments);
    }
    print_eds(msa,segments,filename + ".eds.txt");
    cout << "Minimum segmentation cardinality after gap removal: " << card_eds(msa, segments) << "\n";
    return 0;
}