#include "rmqueue.h"

#include <numeric>
#include <math.h>
#include <iostream>

using std::cout;

// Simple function for floor of base 2 logarithm
long long log2(long long x) {
  int y = 0; 
  while (x >>= 1) {
    y++;
  } 
  return y;
}

// Ceil division operation
long long ceil_div(long long x, long long y) {
  if (x % y == 0) {
    return x / y;
  } else {
    return x / y + 1;
  }
}

// API

RMQueue::RMQueue(long long n): n(n), logn(log2(n)) {
  b = ceil_div(logn, 4);
  blocks = ceil_div(n, b) + 1;
  q_size = b * blocks;
  Q = vector(q_size, std::numeric_limits<long long>::max());
  // Initialize logarithm array
  lsize = 1 << ceil_div(logn, 2); // O(sqrt(n))
  L.resize(lsize + 1);
  L[0] = 0;
  long long i = 1, lg = 0, count = 1;
  while(i <= lsize) {
    L[i++] = lg;
    count--;
    if(count == 0){
      lg++;
      count = 1 << lg;
    }
  }
  // Block aligned queries
  p = vector<vector<long long>>(logn + 1, vector<long long>(blocks + 1, -1));
  // Precompute tables
  compute_tj();
  compute_lca();
  // In block queries
  cart = vector<long long>(blocks, 0);
  spine = vector<long long>(b + 1, 0);
}

void RMQueue::pop() {
  start++;
  if (start >= b) {
    start = 0;
    begin = plus(begin, 1);
  }
}

void RMQueue::push(long long x) {
  setq(end * b + stop, x);
  // Update C array
  update_cartesian_tree(x);
  stop++;
  if (stop >= b) {
    stop = 0;
    // Update arrays with new block
    update_pk(end);
    // Next block
    end = plus(end, 1);
    // Clear cartesian tree at that block
    cart[end] = 0;
    spine_i = 0;
  }
}

// Get element at index i in Q
long long RMQueue::get(long long i) {
  if (i < 0 || i > n) {
    return -1;
  }
  return getq(i + start + begin * b);
}

long long RMQueue::query(long long left, long long right) {
  if (left < 0 || right > n || left > right) {
    return -1;
  }
  long long i = left + start;
  long long j = right + start;
  // Inclusive with first/last block element
  long long b1 = i / b; // block of i
  long long b2 = j / b; // block of j
  long long offset = begin * b + start;
  long long x = -1;
  if (b1 == b2) {
    x = in_block_query(b1, i % b, j % b);
  } else {
    long long q1 = in_block_query(b1, i % b, b - 1);
    long long q2 = block_aligned_query(b1 + 1, b2 - 1);
    long long q3 = in_block_query(b2, 0, j % b);
    x = minx(minx(q1, q2), q3);
  }
  return (x + q_size - offset) % q_size;
}

// Util methods

long long RMQueue::block_aligned_query(long long b1, long long b2) {
  // Check if query is empty
  if (b2 < b1)
    return -1;
  // Block query [b1, b2]
  long long k = log(b2 - b1);
  long long x1 = p[k][plus(b2, begin)];
  long long x2 = p[k][plus(b1 - 1 + (1 << k), begin)];
  return minx(x2, x1);
}

long long RMQueue::in_block_query(long long block, long long i1, long long i2) {
  if (i1 == i2)
    return plus(block, begin) * b + i1;
  if (plus(block, begin) != end) {
    // Complete block
    return plus(block, begin) * b + lca[cart[plus(block, begin)]][1 + i1][1 + i2] - 1;
  }
  // Last block is partial => cartesian tree with #stop number of vertices
  // Get a valid cartesian tree that starts the same
  long long t = (cart[plus(block, begin)] << ((b - stop) * 2)) + (1 << (b - stop)) - 1;
  return plus(block, begin) * b + lca[t][1 + i1][1 + i2] - 1;
}

long long RMQueue::getq(long long i) {
  return Q[i % q_size];
}

void RMQueue::setq(long long i, long long x) {
  Q[i % q_size] = x;
}

long long RMQueue::minx(long long x1, long long x2) {
  if (x1 == -1)
    return x2;
  if (x2 == -1)
    return x1;
  if(getq(x1) <= getq(x2)) {
    return x1;
  }
  return x2;
}

long long RMQueue::plus(long long x, long long y) {
  return (x + y + blocks * 10) % blocks;
}

long long RMQueue::log(long long x) {
  if(x <= lsize)
    return L[x];
  return L[x / lsize] + ceil_div(logn, 2);
}

void RMQueue::update_pk(long long block) {
  long long min_b = std::numeric_limits<long long>::max();
  long long l = -1; // index
  for (int i = 0; i < b; i++) {
    long long val = getq(block * b + i);
    if (min_b > val) {
      min_b = val;
      l = block * b + i;
    }
    p[0][block] = l;
    for (int k = 1; k <= logn; k++) {
      long long l1 = p[k - 1][block];
      long long l2 = p[k - 1][plus(block, - (1 << (k - 1)))];
      p[k][block] = minx(l2, l1);
    } 
  }
}

void RMQueue::compute_tj() {
  tj.resize(b); // [1,b-1]
  tj[0] = vector<vector<long long>>(1, vector<long long>(1, 1));
  vector<long long> right_spine(b * 2, -1);
  right_spine[0] = 0;
  for (int j = 1; j < b; j++) {
    tj[j].resize(1 << (2 * j));
    // All cartesian trees x with j vertices
    // Balanced parentheses 0 = (, 1 = )
    for (int x = 0; x < (1 << (2 * j)); x++) {
      long long i = 0; // right_spine index
      long long open = 0;
      long long v = 0; // tree node
      bool ok = true;
      for (int bit = 2 * j - 1; bit >= 0; bit--) {
        if (open == -1) {
          ok = false;
          break;
        }
        if (open == 0) {
          // Start new spine
          i = 0;
        }
        if ((x >> bit) & 1) {
          open--;
        } else {
          open++;
          // Add v to spine
          right_spine[++i] = ++v;
        }
      }
      if (ok && open == 0) {
        // Integer x describes a cartesian tree
        tj[j][x] = vector(j + 1, -1LL);
        long long new_x = (x << 2) + 3; // add 2 ) at the end
        for (int i2 = 0; i2 <= i; i2 ++) {
          // set ( at position i2 from the end
          tj[j][x][right_spine[i2]] = new_x & ~(1 << (i2 + 1));
        }
      }
    }
  }
}

int calculate_level(vector<int>& parent, int node) {
  int level = 0;
  while(parent[node] != 0) {
    level++;
    node = parent[node];
  }
  return level;
}

void RMQueue::compute_lca() {
  cout << "Computing LCA b = " << b << "\n";
  lca.resize(1 << (2 * b));
  vector<int> parent(b + 1, -1);
  vector<int> level(b + 1, -1);
  vector<int> stack_trace(b + 1, -1);
  // All cartesian trees with vertices
  for (int x = 0; x < (1 << (2 * b)); x++) {
      long long open = 0;
      long long i = -1; // index for stack-trace
      long long v = 0; // current number of (
      long long last = -1; // last closed node
      bool ok = true;
      for (int bit = 2 * b - 1; bit >= 0; bit--) {
        if ((x >> bit) & 1) {
          open--;
        } else {
          open++;
        }
        if (open < 0 || open > b) {
          ok = false;
          break;
        }
        if ((x >> bit) & 1) {
          // Close last node from stack trace
          if (last != -1) {
            parent[last] = stack_trace[i];
          }
          last = stack_trace[i];
          i--;
        } else {
          // Add next v to trace
          stack_trace[++i] = ++v;
          if (last != -1) {
            parent[last] = v;
            last = -1;
          }
        }
      }
      if (ok && open == 0) {
        // Integer x describes a cartesian tree
        lca[x] = vector(b + 1, vector(b + 1, -1LL));
        // Simply calculate lca of each 2 nodes
        parent[last] = 0; // last node closed is root
        for (int i = 1; i <= b; i++) {
          level[i] = calculate_level(parent, i);
        }
        for (int p = 1; p <= b; p++) {
          for (int q = 1; q <= b; q++) {
            // LCA for p & q
            int parent_p = p;
            int parent_q = q;
            while (parent_p != parent_q) {
              if (level[parent_p] > level[parent_q])
               parent_p = parent[parent_p];
              else 
               parent_q = parent[parent_q];
            }
            lca[x][p][q] = parent_p;
          }
        }
      }
    }
}

void RMQueue::update_cartesian_tree(long long x) {
  // Find where to insert x on the spine
  while(spine_i > 0 && getq(end * b + spine[spine_i] - 1) > x) {
    spine_i--;
  }
  // Update the cartesian tree by inserting the new node
  cart[end] = tj[stop][cart[end]][spine[spine_i]];
  // Place the new node on the right spine
  spine_i++;
  spine[spine_i] = stop + 1;
}

// Debugging

void print_n_bits(long long x, int n) {
  for (int i = n - 1; i >= 0; i--) {
    cout << ((x >> i) & 1LL);
  }
}

void RMQueue::debug_tables() {
  cout << "Debugging Tj table:\n";
  for (int j = 1; j < b; j++) {
    cout << "Tj[" << j << "]:\n";
    for (int x = 0; x < (1 << (2 * j)); x++) {
      if(!tj[j][x].empty()) {
        cout << "x = ";
        print_n_bits(x, j * 2);
        cout << "\n";
        for (int h = 0; h <= j; h++) {
          if (tj[j][x][h] != -1) {
            cout << "  " << h << " -> ";
            print_n_bits(tj[j][x][h], j * 2 + 2);
            cout << "\n";
          }
        }
      }
    }
  }
  cout << "Debugging LCA table:\n";
  for (int x = 0; x < (1 << (2 * b)); x++) {
    if(!lca[x].empty()) {
      cout << "LCA for x = ";
      print_n_bits(x, b * 2);
      cout << "\n";
      for (int p = 1; p <= b; p++) {
        for (int q = p; q <= b; q++) {
          cout << p << " " << q << " -> " << lca[x][p][q] << "\n";
        }
      }
    }
  }
}

void RMQueue::debug_trees() {
  cout << "Debugging cartesian trees:\n";
  if (end >= begin) {
    // No ring buffer wrapping
    for (int i = begin; i < end; i++) {
      cout << "c[" << i << "] = ";
      print_n_bits(cart[i], 2 * b);
      cout << "\n";
    }
  } else {
    for (int i = begin; i < blocks; i++) {
      cout << "c[" << i << "] = ";
      print_n_bits(cart[i], 2 * b);
      cout << "\n";
    }
    for (int i = 0; i < end; i++) {
      cout << "c[" << i << "] = ";
      print_n_bits(cart[i], 2 * b);
      cout << "\n";
    }
  }
  // Current block
  cout << "current cartesian tree ";
  print_n_bits(cart[end], 2 * stop);
  cout << "\n";
}