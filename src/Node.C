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
	this->occs = tmpNode.occs.load();
	this->weight = tmpNode.weight.load();
	this->endCnt = tmpNode.endCnt.load();
	omp_init_lock( &lock );
	this->offset = tmpNode.offset;
	this->readNum = tmpNode.readNum;

	for( int i = 0; i < NBASES; i++ )
		this->subnodes[i] = tmpNode.subnodes[i];
}

Node& Node::operator=( const Node& tmpNode ) 
{
	this->occs = tmpNode.occs.load();
	this->weight = tmpNode.weight.load();
	this->endCnt = tmpNode.endCnt.load();
	omp_init_lock( &lock );
	this->offset = tmpNode.offset;
	this->readNum = tmpNode.readNum;

	for( int i = 0; i < NBASES; i++ )
		this->subnodes[i] = tmpNode.subnodes[i];

	return* this;
}

double Node::getRatio()
{
	return weight.load() / static_cast<double>( occs.load() );
}
