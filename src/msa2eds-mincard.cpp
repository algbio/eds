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
#include <tuple>
#include <iomanip>

#include "RMaxQTree.h"
#include "block_graph.hpp"
#include "msa_chunker.hpp"
#include "CLI11.hpp"

using namespace std::chrono;
using namespace std;
using eds::block_graph::block_graph, eds::block_graph::segment_msa, eds::block_graph::output_msa_info, eds::block_graph::output_segmentation, eds::block_graph::output_block_info, eds::block_graph::output_block_graph, eds::block_graph::output_eds;
using msa_chunker::msa_chunker, msa_chunker::msa_pos_t;

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

vector<pair<seg_index, seg_index>> compute_meaningful_extensions(
    msa_chunker::msa_chunker &idx,
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

    // Enforce ℓ_{y,1} = y - L + 1 down to ℓ_{y,d_y} > y - U
    seg_index prev_height = -1;
    for (seg_index len = L; len <= U && y - len + 1 >= 1; ++len) {
        seg_index start = y - len + 1;
        unordered_set<string> unique_strings;

        for (seg_index i = 0; i < r; ++i) {
            string s = idx.msa_substr(i, start - 1, len);
            if (!gaps_as_symbols)
              s.erase(remove(s.begin(), s.end(), '-'), s.end());
            unique_strings.insert(s);
        }

        seg_index height = unique_strings.size();
        if (height != prev_height) {
            L_y.emplace_back(start, height);
            prev_height = height;
        }
    }

    // Add dummy ℓ_{y,d_y+1} = max(0, y - U)
    L_y.emplace_back(max((seg_index)0, y - U), -1);

    return L_y;
}
vector<vector<pair<seg_index, seg_index>>> compute_all_meaningful_extensions(
    msa_chunker::msa_chunker &idx,
    const seg_index r,
    const seg_index c,
    const seg_index L,
    const seg_index U,
    const bool gaps_as_symbols
) {
    vector<vector<pair<seg_index, seg_index>>> L_y(c + 1);  // 1-based indexing

    for (seg_index y = 1; y <= c; ++y) {
        L_y[y] = compute_meaningful_extensions(idx, r, c, L, U, y, gaps_as_symbols);
    }

    return L_y;
}

pair<seg_index,vector<bool>> compute_perfect_columns(
    msa_chunker::msa_chunker &idx, const seg_index r, const seg_index c) {
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
pair<seg_index, vector<pair<seg_index, seg_index>>> segment_with_rmq(
        msa_chunker::msa_chunker &idx,
        const seg_index r,
        const seg_index c,
        const int L,
        const int U,
        const bool gaps_as_symbols,
        const vector<vector<pair<seg_index, seg_index>>> &L_y,
        const vector<bool> &perfect_columns = perfect_columns_dummy
) {
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

        //const auto& L = L_y[y];
        const vector<pair<seg_index, seg_index>> *L_yy;
        vector<pair<seg_index, seg_index>> L_yy_on_the_fly;
        if (L_y.size() > 0) {
          L_yy = &(*(L_y.begin() + y));
        } else {
          L_yy_on_the_fly = compute_meaningful_extensions(idx, r, c, L, U, y, gaps_as_symbols);
          L_yy = &L_yy_on_the_fly;
        }

        // optimal solution using L_yy
        for (size_t j = 0; j + 1 < L_yy->size(); ++j) {
            key_type l = (*L_yy)[j + 1].first;
            key_type r = (*L_yy)[j].first - 1;
            if (l > r) continue;

            // query returns pair (index, value), but value is -m[index]
            auto [x, neg_mx] = rmq.query(l, r);
            key_type candidate = (*L_yy)[j].second + m[x];

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
void prseg_index_eds(msa_chunker::msa_chunker &idx, const vector<pair<seg_index, seg_index>>& segments, string out_filename = "") {
    std::ofstream outFile;
    const seg_index rows = idx.get_rows();
    if (out_filename.size()==0)
        cout << "Elastic Degenerate String (EDS):\n";
    else
       outFile.open(out_filename);
    for (const auto& [l, r] : segments) {
        set<string> unique_subs;
        for (seg_index i = 0; i < rows; ++i) {
            string sub = idx.msa_substr(i, l - 1, r - l + 1);
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
seg_index card_eds(msa_chunker::msa_chunker &idx, const vector<pair<seg_index, seg_index>>& segments) {
    const seg_index rows = idx.get_rows();
    seg_index card = 0;
    for (const auto& [l, r] : segments) {
        set<string> unique_subs;
        for (seg_index i = 0; i < rows; ++i) {
            string sub = idx.msa_substr(i, l - 1, r - l + 1);
            sub.erase(remove(sub.begin(), sub.end(), '-'), sub.end()); // remove gaps
            unique_subs.insert(sub);
        }
        card += unique_subs.size();
    }
    return card;
}

int main(int argc, char* argv[]) {
    CLI::App app{"msa2eds-mincard version " + string(VERSION) + " — build Elastic Degenerate Strings (EDSes) from multiple sequence alignments in FASTA format"};
    argv = app.ensure_utf8(argv);

    string inputfile;
    app.add_option("msa.fasta", inputfile, "Input MSA (FASTA format)")
      ->required();

    string outputedsfile;
    app.add_option("-o,--output-eds", outputedsfile, "Output file for the elastic degenerate string (EDS format)")
      ->default_val("");

    string outputgfafile;
    app.add_option("-g,--output-gfa", outputgfafile, "Output file for the corresponding block graph (xGFA format, not recommended with --perfect-segments)")
      ->default_val("");

    seg_index L;
    CLI::Option *Lopt = app.add_option("-L,--min-segment-length", L, "Minimum segment length")
      ->default_val(1)
      ->expected(1, msa_chunker::MAX_CHUNK_COLS);

    seg_index U;
    CLI::Option *Uopt = app.add_option("-U,--max-segment-length", U, "Maximum segment length")
      ->default_val(31)
      ->expected(1, msa_chunker::MAX_CHUNK_COLS);

    bool allow_perfect_segments = false;
    app.add_flag("-p,--perfect-segments", allow_perfect_segments, "In normal mode, additionally consider perfect segments of any length (recommended). With --trivial-vertical and --trivial-horizontal, use the maximal perfect segments and the trivial strategy in-between.");

    bool trivial_segmentation = false;
    CLI::Option *tsopt = app.add_flag("-t,--trivial-vertical", trivial_segmentation, "Use trivial S^¦¦¦ segmentation (every column becomes an ED word)")
      ->excludes(Lopt)->excludes(Uopt);

    bool no_segmentation = false;
    CLI::Option *nsopt = app.add_flag("-n,--trivial-horizontal", no_segmentation, "Use trivial S^≡ segmentation (no segmentation)")
      ->excludes(Lopt)->excludes(Uopt)->excludes(tsopt);

    bool gaps_as_symbols = false;
    app.add_flag("--gaps-as-symbols", gaps_as_symbols, "In preprocessing the MSA, consider gaps '-' as normal symbols")
      ->excludes(tsopt)->excludes(nsopt);

    bool preprocess = false;
    app.add_flag("--preprocess", preprocess, "Compute all meaningful extensions before segmenting")
      ->excludes(tsopt)->excludes(nsopt);

    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "Print running times");

    try {
      app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
      return app.exit(e);
    }

    msa_chunker::msa_chunker idx(inputfile, U);
    const int r = idx.get_rows();
    const int c = idx.get_cols();
    cerr << "Processing MSA[1.." << r << ",1.." << c << "] (\"" << inputfile << "\")" << endl;

    vector<bool> perfect_columns = {};
    if (allow_perfect_segments) {
      seg_index p;
      auto start = high_resolution_clock::now();
      tie(p, perfect_columns) = compute_perfect_columns(idx, r, c);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cout << "MSA contains " << p << "/" << c << " (" << setprecision(4) << (double) 100 * p / c << "%) perfect columns" <<((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    }

    vector<pair<seg_index, seg_index>> segmentation; // 1-based segments [x..y]
    if (trivial_segmentation) {
      cerr << "Computing the S^¦¦¦ segmentation..." << flush;
      segmentation.reserve(c);
      for (seg_index i = 1; i <= c; ++i) {
        seg_index j = i;
        if (allow_perfect_segments and perfect_columns[j]) {
          while (perfect_columns[j]) 
            j += 1;
          j -= 1;
        }
        segmentation.push_back({ i, j });
        i = j;
      }
      cerr << " done: "  << segmentation.size() << " segments/ED words" << endl;
    } else if (no_segmentation) {
      cerr << "Computing the S^≡ segmentation..." << flush;
      if (!allow_perfect_segments) {
        segmentation.push_back({ 1, c });
      } else {
        for (seg_index i = 1; i <= c; ++i) {
          seg_index j = i;
          bool start = perfect_columns[j];
          while (j <= c and perfect_columns[j] == start) {
            j += 1;
          }
          j -= 1;
          segmentation.push_back({ i, j });
          i = j;
        }
      }
      cerr << " done: "  << segmentation.size() << " segments/ED words" << endl;
    } else {
      cerr << "The allowed segments are" << ((allow_perfect_segments) ? " perfect segments and those" : "") << " of length [" << L << ".." << U << "]" << endl;

      vector<vector<pair<seg_index, seg_index>>> L_y;
      if (preprocess) {
	      cerr << "Computing the meaningful extensions..." << flush;
	      auto start = high_resolution_clock::now();
	      L_y = compute_all_meaningful_extensions(idx, r, c, L, U, gaps_as_symbols);
	      auto stop = high_resolution_clock::now();
	      auto duration = duration_cast<milliseconds>(stop - start);
	      cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
      }

      cerr << "Computing the minimum-cardinality segmentation..." << flush;
      auto start = high_resolution_clock::now();
      seg_index mincard;
      tie(mincard, segmentation) = segment_with_rmq(idx, r, c, L, U, gaps_as_symbols, L_y, perfect_columns);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cerr << " done: " << segmentation.size() << " segments/ED words, " << mincard << " cardinality" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    }

    eds::block_graph::seg_size_t card, size;
    if (outputgfafile != "") {
      ostream *out;
      ofstream outfile;
      if (outputgfafile == "-") {
        cerr << "Streaming the block graph to stdout..." << flush;
        out = &cout;
      } else {
        cerr << "Streaming the block graph to \"" << outputgfafile << "\"..." << flush;
        outfile = ofstream(outputgfafile);
        out = &outfile;
      }
      auto start = high_resolution_clock::now();
      tie(card, size) = eds::block_graph::segment_stream_gfa(idx, r, c, segmentation, out);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      outfile.close();
      cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << ":";
    }
    if (outputedsfile != "") {
      ostream *out;
      ofstream outfile;
      if (outputedsfile == "-") {
        cerr << "Streaming the EDS to stdout..." << flush;
        out = &cout;
      } else {
        cerr << "Streaming the EDS to \"" << outputedsfile << "\"..." << flush;
        outfile = ofstream(outputedsfile);
        out = &outfile;
      }
      auto start = high_resolution_clock::now();
      tie(card, size) = eds::block_graph::segment_stream_eds(idx, r, c, segmentation, out);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      outfile.close();
      cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << ":";
    }
    if (outputgfafile == "" and outputedsfile == "") {
      cerr << "Computing the EDS stats (no output selected)..." << flush;
      auto start = high_resolution_clock::now();
      tie(card, size) = eds::block_graph::segment_stream_no_output(idx, r, c, segmentation);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cerr << " done" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << ":";
    }

    cerr << " " << card << " cardinality, " << size << " gap-aware size" << endl;
    return 0;
}
