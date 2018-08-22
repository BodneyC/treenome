/********************************************************************
 * Filename: Node.H [C++ header code]
 *
 * Description: Header for Node class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
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
	uint64_t readNum;
	int32_t subnodes[NBASES];

	Node();
	//Copy and assignment functions needed for atomics
	Node( const Node& tmpNode );
	Node& operator=( const Node& tmpNode );

	double getRatio();
} Node;

#endif /*_NODE*/
