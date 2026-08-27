#include "pbwt.h"

#include <numeric>
#include <optional>
#include <algorithm>
#include "rmq.hpp"

pbwt::pbwt(const long long r): 
  r(r), k(0),
  ak(r + 1, 0), sk(r + 2, 0), tk(r + 2, 0), ek(r + 1, 0), 
  dy(0), a(r + 2, 0), e(r + 1, 0)
{
  cnt.resize(alphabet_size + 1);
  prev.resize(alphabet_size);
  // Initial sorted order 1, 2, 3..
  std::iota(ak.begin(), ak.end(), 0);
}

// If the symbol doesn't exist in the alpabet, add it
unsigned int pbwt::symbol_number(const char symbol) {
  auto s = static_cast<unsigned char>(symbol);

  if (alphabet[s] == -1) {
    alphabet[s] = alphabet_size++;
    cnt.push_back(0);
    prev.push_back(0);
  }

  return alphabet[s];
}

long long pbwt::max_ek(long long j, long long i) {
  if (j != i) {
    ek[j] = std::max(ek[j], max_ek(ak[j], i));
    ak[j] = i + 1;
  }
  return ek[j];
}

void pbwt::next(const string &column) {
  // Increment column
  k++;
  // Symbol frequency
  cnt.assign(cnt.size(), 0);
  for (long long i = 0; i < r; i++){
    auto c = column[i];
    auto s = symbol_number(c) + 1;
    cnt[s]++;
  }
  int symbols = 0;
  for (long long sym = 1; sym <= alphabet_size; sym++)
    if(cnt[sym] > 0)
      symbols++;
  if (k == 1 or symbols > 1) {
    // There is at least one new divergence => height change
    // Recompute pBWT arrays
    prev.assign(prev.size(), 0);
    tk.assign(tk.size(), 0);
    sk[++dy] = k;
    // RMQ for max ek
    std::optional<rmq<long long>> rm;
    if constexpr (max_range == RMQ) {
      rm.emplace(ek);
    }
    // Counting sort
    for (unsigned int i = 1; i < alphabet_size; i++)
      cnt[i] += cnt[i - 1];
    for (long long i = 1; i <= r; i++) {
      unsigned int b = symbol_number(column[ak[i] - 1]);
      cnt[b]++;
      a[cnt[b]] = ak[i];
      if constexpr (max_range == RECURSIVE) {
        ak[i] = i + 1;
      }
      if (prev[b] == 0) 
        e[cnt[b]] = dy;
      else {
        // Calculate the range maximum 
        if constexpr (max_range == NAIVE) {
          // O(alphabet) amortized - very fast in practice
          e[cnt[b]] = *std::max_element(ek.begin() + prev[b] + 1, ek.begin() + i + 1);
        } else if constexpr (max_range == RECURSIVE) {
          // O(log alphabet) - solution from the paper
          e[cnt[b]] = max_ek(prev[b] + 1, i);
        } else {
          // RMQ O(1) - best complexity but longer construction
          e[cnt[b]] = rm->query(prev[b] + 1, i);
        }
      }
      prev[b] = i;
    }
    for (long long i = 1; i <= r; i++) {
      ak[i] = a[i];
    }
    // Calculate frequency array
    for (long long i = 1; i <= r; i++) {
      tk[e[i]]++;
    }
    // Shrink arrays sk and tk, array a acts as tmp
    long long j = 1;
    for (long long i = 1; i <= dy; i++) {
      if (tk[i] != 0) {
        a[i] = j;
        sk[j] = sk[i];
        tk[j] = tk[i];
        j++;
      }
    }
    dy = j - 1;
    // Fix ek array
    for (long long i = 1; i <= r; i++) {
      ek[i] = a[e[i]];
    }
    // Calculate heights of extensions
    for(long long i = 1; i <= dy; i++){
      tk[dy - i] += tk[dy - i + 1];
    }
  }
  else {
    // Update last element of sk for first row divergence
    if (tk[dy] > 1){
      dy++;
    }
    sk[dy] = k;
    ek[1] = dy;
    tk[dy] = 1;
  }
}

vector<pair<long long, long long>> pbwt::meaningful_extensions(
  const long long L,
  const long long U
) {
  if (k < L) {
    // No extension possible
    return {};
  }
  vector<pair<long long, long long>> L_y; 
  for(long long i = 0; i < dy; i++){
    if (k - sk[dy - i] + 1 < L) {
      if (k - sk[dy - i - 1] + 1 > L) {
        // Introduce new left extension at k - L
        L_y.emplace_back(k - L + 1, tk[dy - i]);       
      }
    } else {
      if (k - sk[dy - i] + 1 <= U) {
        // Left extension is in range [k - U, k - L]
        L_y.emplace_back(sk[dy - i], tk[dy - i]);
      }
    }
  }
  // Dummy element, extension at k - U - 1
  L_y.emplace_back(std::max((long long)0, k - U), -1);
  return L_y;
}

long long pbwt::get_column(){
  return k;
}
vector<long long> pbwt::get_ak(){
  return ak;
}
vector<long long> pbwt::get_dk(){
  vector<long long> d(r + 1, 0);
  for(int i = 1; i <= r; i++) {
    d[i] = sk[ek[i]];
  }
  return d;
}
