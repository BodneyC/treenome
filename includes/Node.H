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

typedef struct BNode {
	std::atomic<int64_t> occs;
	std::atomic<double> weight;
	std::atomic<int32_t> endCnt;
	omp_lock_t lock;

	BNode();
	//Copy and assignment functions needed for atomics
	BNode( const BNode& tmpNode );
	BNode& operator=( const BNode& tmpNode );

	double getRatio();
} BNode;

typedef struct Node : public BNode {
	short offset;
	int64_t readNum;
	Node* subnodes[NBASES];

	Node();
} Node;

typedef struct pNode : BNode {
	pNode* subnodes[NBASES];
	pNode* parent;

	pNode();
} pNode;

#endif /*_NODE*/
