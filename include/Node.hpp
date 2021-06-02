/********************************************************************
 * Filename: Node.hpp [C++ header code]
 *
 * Description: Header for Node class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#ifndef _NODE_
#define _NODE_

#if __cplusplus < 201103L
#error This code requires C++11
#endif /*c++11 check*/

#define NBASES 4

#include <atomic>
#include <omp.h>

extern int NUM_THREADS;

typedef struct Node {
  std::atomic<int64_t> occs;
  std::atomic<double> weight;
  omp_lock_t lock;
  short offset;
  uint64_t read_num;
  int32_t subnodes[NBASES];

  Node();
  // Copy and assignment functions needed for atomics
  Node(const Node &tmpNode);
  Node &operator=(const Node &tmpNode);

  double get_ratio();
} Node;

#endif /*_NODE*/
