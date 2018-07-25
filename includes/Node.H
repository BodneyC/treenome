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
#error This code requires C++11 for threading
#endif /*c++11 check*/

#define NBASES 4
#define RES 50
#define NUM_THREADS 4
#include <atomic>
#include <omp.h>

// The node type for the GTrees 
typedef struct Node {
	int64_t occs;
	std::atomic<double> weight;
	short offset;
	int64_t readNum;
	Node* subnodes[NBASES];
	omp_lock_t lock;

	Node();
	//Copy and assignment functions needed for atomics
	Node(const Node& tmpNode);
	Node& operator=(const Node& tmpNode);
} Node;

typedef struct pNode : Node {
	Node* parent;

	pNode(): parent(nullptr) {  }
} pNode;

#endif /*_NODE*/
