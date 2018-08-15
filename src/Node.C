/********************************************************************
 * Filename: Node.C [C++ source code]
 *
 * Description: Implementation of Node class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/Node.H"

Node::Node(): occs( 0 ), weight( 0 ), endCnt( 0 )
{
	omp_init_lock( &lock );
	for( int i = 0; i < NBASES; i++ )
		subnodes[i] = 0;
}

Node::Node( const Node& tmpNode ) 
{
	float tmpWeight = tmpNode.weight;
	int64_t tmpOccs = tmpNode.occs;

	omp_init_lock( &lock );
	weight = tmpWeight;
	occs = tmpOccs;
}

Node& Node::operator=( const Node& tmpNode ) 
{
	float tmpWeight = tmpNode.weight;
	int64_t tmpOccs = tmpNode.occs;

	omp_init_lock( &lock );
	this->weight = tmpWeight;
	this->occs = tmpOccs;

	return* this;
}

double Node::getRatio()
{
	return weight / static_cast<double>( occs );
}

