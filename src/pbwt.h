#ifndef PBWT_H
#define PBWT_H

#include <vector>
#include <array>
#include <iostream>

using std::vector;
using std::string;
using std::pair;
using std::array;
using std::cerr, std::endl;

enum MaxRange {
    NAIVE,
    RECURSIVE,
    RMQ
};

inline constexpr enum MaxRange max_range = NAIVE;

class pbwt {
  private:
    // Number of rows
    long long r;
    // Current column
    long long k = 0;
    // Arrays for current column 1-based indexing
    vector<long long> ak;
    vector<long long> sk;
    vector<long long> tk;
    vector<long long> ek;
    // Counting sort arrays
    vector<long long> cnt;
    vector<long long> prev;
    // Extra space
    long long dy; // size of sk
    vector<long long> a;
    vector<long long> e;
    // Alphabet (may grow with new symbols)
    array<int, 256> alphabet = [] {
      array<int, 256> alph{};
      alph.fill(-1);

      alph['A'] = 0;
      alph['C'] = 1;
      alph['G'] = 2;
      alph['T'] = 3;
      alph['-'] = 4;

      return alph;
    }();
    unsigned int alphabet_size = 5;
    // Map char -> int
    unsigned int symbol_number(const char symbol);
    // Recursive function for range max
    long long max_ek(long long j, long long i);
  public:
    // Constructor taking the number of rows
    pbwt(const long long r);
    ~pbwt() = default;
    // 1 indexed algorithm components
    long long get_column();
    vector<long long> get_ak();
    vector<long long> get_dk();
    // Column streaming algorithm
    void next(const string &column);
    vector<pair<long long, long long>> meaningful_extensions(
      const long long L,
      const long long U
    );
};

#endif
