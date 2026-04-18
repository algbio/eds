#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <limits>
#include <tuple>
#include <iomanip>

#include "algo.hpp"
#include "segment.hpp"
#include "msa_chunker.hpp"
#include "CLI11.hpp"

using std::cout, std::cerr;
using std::vector;
using std::string, std::to_string;
using std::chrono::high_resolution_clock, std::chrono::duration_cast, std::chrono::milliseconds;
using std::numeric_limits;
using std::pair, std::tie;
using std::setprecision;
using algo::seg_index, algo::compute_perfect_columns, algo::compute_all_meaningful_extensions, algo::segment_with_rmq;
typedef segment::seg_size_t seg_size_t;
using segment::segment_stream_gfa, segment::segment_stream_eds, segment::segment_stream_no_output;

bool verbose = false;

int main(int argc, char* argv[]) {
    CLI::App app{"mincard version " + string(VERSION) + " — build Elastic Degenerate Strings (EDSes) from multiple sequence alignments in FASTA format"};
    argv = app.ensure_utf8(argv);

    string inputfile;
    app.add_option("msa.fasta", inputfile, "Input MSA (FASTA format)")
      ->required();

    string outputedsfile;
    app.add_option("-o,--output-eds", outputedsfile, "Output file for the elastic degenerate string (EDS format)")
      ->default_val("");

    string outputgfafile;
    app.add_option("-g,--output-gfa", outputgfafile, "Output file for the corresponding block graph (xGFA format, not recommended with --trivial-vertical)")
      ->default_val("");

    seg_index L;
    CLI::Option *Lopt = app.add_option("-L,--min-segment-length", L, "Minimum segment length")
      ->default_val(1)
      ->expected(1, numeric_limits<int>::max());

    seg_index U;
    CLI::Option *Uopt = app.add_option("-U,--max-segment-length", U, "Maximum segment length")
      ->default_val(31)
      ->expected(1, numeric_limits<int>::max());

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
    if (L > U) {
      cerr << "Upper and lower bounds are not compatible!" << endl;
      return 1;
    }

    msa_chunker::msa_chunker idx(inputfile, U, verbose);
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
      cerr << "MSA contains " << p << "/" << c << " (" << setprecision(4) << (double) 100 * p / c << "%) perfect columns" <<((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    }

    vector<pair<seg_index, seg_index>> segmentation; // 1-based segments [x..y]
    if (trivial_segmentation) {
      cerr << "Computing the S^¦¦¦ segmentation..." << flush;
      auto start = high_resolution_clock::now();
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
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cerr << " done: "  << segmentation.size() << " segments/ED words" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    } else if (no_segmentation) {
      cerr << "Computing the S^≡ segmentation..." << flush;
      auto start = high_resolution_clock::now();
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
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cerr << " done: "  << segmentation.size() << " segments/ED words" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
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
      if (mincard != std::numeric_limits<seg_index>::max()) {
        cerr << " done: " << segmentation.size() << " segments/ED words, " << mincard << " cardinality" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
      } else {
        cerr << " done: no valid segmentation found!" << endl;
        return 1;
      }
    }

    seg_size_t card, size;
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
      tie(card, size) = segment_stream_gfa(idx, r, c, segmentation, out);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      outfile.close();
      cerr << " done: " << card << " cardinality, " << size << " gap-aware size" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
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
      tie(card, size) = segment_stream_eds(idx, r, c, segmentation, out);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      outfile.close();
      cerr << " done: " << card << " cardinality, " << size << " gap-aware size" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    }
    if (outputgfafile == "" and outputedsfile == "") {
      cerr << "Computing the EDS stats (no output selected)..." << flush;
      auto start = high_resolution_clock::now();
      tie(card, size) = segment_stream_no_output(idx, r, c, segmentation);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      cerr << " done: " << card << " cardinality, " << size << " gap-aware size" << ((verbose) ? " (" + to_string(duration.count()) + "ms)" : "") << endl;
    }

    return 0;
}
