/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree subclass, for building the 
 *		tree from a file outputted from the program
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/GTree.H"

void GTree::getNextNode( struct NodeInfo& nInf, std::stringstream& ss )
{
	std::string line;
	char tmpChar;
	while( 1 ) {
		if( !(ss >> tmpChar) )
			return;
		if( tmpChar != ',' )
			break;
		nInf.comCnt++;
	}
	nInf.ind = BASE_IND( tmpChar );
	ss >> tmpChar;
	std::getline( ss, line, ':' );
	nInf.occs = std::stol( line );
	std::getline( ss, line, ';' );
	nInf.weight = std::stod( line );
}

void GTree::createNode( struct NodeInfo& nInf )
{
	for( int i = 0; i < nInf.comCnt; i++ ) {
		//std::cout << tmpNode << std::endl;
		tmpNode = tmpNode->parent;
	}
	
	tmpNode->subnodes[nInf.ind] = &( vNodes[head] );
	head++;

	tmpNode->subnodes[nInf.ind]->parent = tmpNode;
	tmpNode = tmpNode->subnodes[nInf.ind];

	tmpNode->occs = nInf.occs;
	tmpNode->weight = nInf.weight;
}

void GTree::createRoot( std::stringstream& ss )
{
	struct NodeInfo nInf;
	getNextNode( nInf, ss );
	root = &( vNodes[0] );
	tmpNode = root;
	head++;

	tmpNode->occs = nInf.occs;
	tmpNode->weight = nInf.weight;
}

void GTree::processSString( std::stringstream& ss )
{
	std::string line;
	std::getline( ss, line, '\n' );

	nNodes = std::stol( line );
	vNodes.resize( nNodes );

	createRoot( ss );

	for( unsigned int i = 0; i < nNodes - 1; i++ ) {
		struct NodeInfo nInf;
		getNextNode( nInf, ss );
		createNode( nInf );
	}
}
