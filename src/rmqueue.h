#ifndef RMQUEUE_H
#define RMQUEUE_H

#include <vector>

using std::vector;

/* Range Minimum Queue
 * mantains a queue of at most N integers and supports operations:
 * pop in O(1)
 *  - deletes the first element of the queue
 * push (element) in amortized O(1)
 *  - inserts the element at the back of the queue
 * query [left, right] in O(1)
 *  - returns the minimum in the specified interval
 */
class RMQueue {
  private:
    // Ring buffer for the actual queue elements
    vector<long long> Q;
    long long getq(long long i);
    void setq(long long i, long long x);
    long long minx(long long x1, long long x2); // index of minimum from 2 ind
    // Invariants
    long long n;
    long long logn;
    long long b; // block length
    long long blocks;
    long long q_size;
    // Current element
    long long begin = 0; // block
    long long start = 0; // shift
    long long end = 0; // block
    long long stop = 0; // in-block shift
    // Modulo addition
    long long plus(long long x, long long y);
    // Precomputed logarithm array
    long long lsize;
    vector<long long> L;
    long long log(long long x);
    // Queries
    long long block_aligned_query(long long b1, long long b2);
    long long in_block_query(long long block, long long i1, long long i2);
    // Block-aligned
    // p[k][h] is the index of the minimum element in Q[(h-2^k)b+1, hb]
    vector<vector<long long>> p; 
    void update_pk(long long block);
    // In-block
    vector<vector<vector<long long>>> tj;
    void compute_tj();
    vector<vector<vector<long long>>> lca;
    void compute_lca();
    vector<long long> cart;
    vector<long long> spine;
    long long spine_i = 0;
    void update_cartesian_tree(long long x);
  public:
    RMQueue(long long n);
    void pop();
    void push(long long x);
    long long query(long long left, long long right);
    long long get(long long i);
    void debug_tables();
    void debug_trees();
};

#endif