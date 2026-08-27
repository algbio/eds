#include "rmqueue.h"
#include <iostream>
#include <assert.h>

using std::cout;

void test_b_1() {
  // Small queue
  RMQueue rmq(7);
  rmq.push(100);
  rmq.push(12);
  rmq.push(50);
  rmq.push(75);
  rmq.push(7);
  // rmq can be indexed
  assert(rmq.get(0) == 100);
  assert(rmq.get(1) == 12);
  assert(rmq.get(4) == 7);
  // test queries
  assert(rmq.query(0, 4) == 4);
  assert(rmq.query(0, 3) == 1);
  assert(rmq.query(2, 2) == 2);
  assert(rmq.query(3, 4) == 4);
  // remove element 100
  rmq.pop();
  assert(rmq.get(0) == 12);
  assert(rmq.query(0, 3) == 3);
  assert(rmq.query(0, 2) == 0);
  assert(rmq.query(1, 2) == 1);
  // add more elements
  rmq.push(21);
  rmq.push(5);
  rmq.pop();
  rmq.push(2014);
  rmq.pop();
  rmq.pop();
  rmq.push(-3);
  rmq.push(9);
  rmq.pop();
  rmq.push(90);
  rmq.push(6);
  // queue 21 5 2014 -3 9 90 6
  assert(rmq.get(0) == 21);
  assert(rmq.get(1) == 5);
  assert(rmq.get(6) == 6);
  assert(rmq.query(0, 6) == 3);
  assert(rmq.query(0, 2) == 1);
  assert(rmq.query(4, 5) == 4);
}

void test_b_3() {
  RMQueue rmq(1<<12);
  // rmq.debug_tables();

  for (int i = 0; i < 600; i++)
    rmq.push(i);
  
  assert(rmq.query(1, 32) == 1);
  assert(rmq.query(20, 400) == 20);
  assert(rmq.query(3, 4) == 3);

  for (int i = 0; i < 300; i++)
    rmq.pop();

  assert(rmq.get(0) == 300);
  assert(rmq.get(rmq.query(2, 5)) == 302);

  for (int i = 600; i < (1<<12); i++)
    rmq.push(i);

  for (int i = 0; i < 200; i++)
    rmq.pop();

  for (int i = 1<<12; i < (1<<12) + 250; i++)
    rmq.push(i); 
    
  assert(rmq.query(0, 12) == 0);
  assert(rmq.query(25, 26) == 25);
  assert(rmq.get(rmq.query(25, 26)) == 525);

  for (int i = 0; i < 3500 ; i++)
    rmq.pop(); 

  assert(rmq.query(14, 122) == 14);
  assert(rmq.get(rmq.query(14, 122)) == 4014);
}

void test_equal_mins() {
  RMQueue rmq(10);
  rmq.push(622);
  rmq.push(620);
  rmq.push(620);
  rmq.push(620);
  rmq.push(621);
  rmq.push(622);
  rmq.push(620);
  assert(rmq.query(0, 0) == 0);
  assert(rmq.query(0, 6) == 1); // first minimum index
}

int main() {
  cout << "Basic RMQueue Testing\n";
  test_b_1();
  test_b_3();
  test_equal_mins();
  cout << "All tests passed!\n";
}