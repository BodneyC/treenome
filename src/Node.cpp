/********************************************************************
 * Filename: Node.cpp [C++ source code]
 *
 * Description: Implementation of Node class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#include "Node.hpp"

Node::Node(): occs( 0 ), weight( 0 )
{
	omp_init_lock( &lock );
	for( int i = 0; i < NBASES; i++ )
		subnodes[i] = 0;
}

Node::Node( const Node& tmpNode ) 
{
	this->occs = tmpNode.occs.load();
	this->weight = tmpNode.weight.load();
	omp_init_lock( &lock );
	this->offset = tmpNode.offset;
	this->read_num = tmpNode.read_num;

	for( int i = 0; i < NBASES; i++ )
		this->subnodes[i] = tmpNode.subnodes[i];
}

Node& Node::operator=( const Node& tmpNode ) 
{
	this->occs = tmpNode.occs.load();
	this->weight = tmpNode.weight.load();
	omp_init_lock( &lock );
	this->offset = tmpNode.offset;
	this->read_num = tmpNode.read_num;

	for( int i = 0; i < NBASES; i++ )
		this->subnodes[i] = tmpNode.subnodes[i];

	return* this;
}

double Node::get_ratio()
{
	return weight.load() / static_cast<double>( occs.load() );
}
