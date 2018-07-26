/********************************************************************
 * Filename: BNode.C [C++ source code]
 *
 * Description: Implementation of BNode class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/Node.H"

BNode::BNode(): occs(0), weight(0)
{
	omp_init_lock(&lock);
}

BNode::BNode(const BNode& tmpBNode) 
{
	float tmpWeight = tmpBNode.weight;
	int64_t tmpOccs = tmpBNode.occs;

	omp_init_lock(&lock);
	weight = tmpWeight;
	occs = tmpOccs;
}

BNode& BNode::operator=(const BNode& tmpBNode) 
{
	float tmpWeight = tmpBNode.weight;
	int64_t tmpOccs = tmpBNode.occs;

	omp_init_lock(&lock);
	this->weight = tmpWeight;
	this->occs = tmpOccs;

	return* this;
}

Node::Node(): offset(0), readNum(0)
{
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = nullptr;
}


pNode::pNode(): parent(nullptr)
{
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = nullptr;
}

