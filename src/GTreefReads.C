/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree subclass, for building the 
 *		tree from a set of reads
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/GTree.H"

void GTree::init() {
	// Because reallocation of std::vectors is assured if its class is a 
	// std::vector, a vector of vectors has been used with a 1-by-1 
	// .push_back()
	omp_init_lock( &lock );
	dNodes.resize( RES );
}
/** ---------------- Tree Creation ----------------- **/
void GTree::updateWeight( Node* node, char qual ) {
	double newWeight, curWeight = node->weight;
	double pBQual = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
	do {
		newWeight = curWeight + pBQual;
	} while( !( node->weight.compare_exchange_weak( curWeight, newWeight ) ) );
}

void GTree::updateWeightAndOccs( Node* node, char qual )
{
	node->occs++;
	updateWeight(node, qual);
}

void GTree::createRoot( short ind )
{
	uint32_t i;
	char qual;
	root = &( dNodes[head] );
	nNodes++;

	for( i = 0; i < GTH::seqReads.size(); i++ ) {
		int offset = -1;
		// Find algorithm in bool vec
		for( int j = 0; j < GTH::seqReads[i].size(); j++ ) {
			if( GTH::seqReads[i].getBaseInd(j) == ind ) {
				offset = j;
				break;
			}
		}
		if( offset != -1 ) {
			qual = GTH::seqReads[i].getQual( offset );
			root->readNum = i;
			root->offset = offset;
			root->weight = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
			root->occs = 1;
			break;
		}
	}
}

void GTree::createNode( Node* node, short ind, char qual, uint64_t rN, int offset )
{
	omp_set_lock( &lock );
	head++;
	if( head == dNodes.size() ) {
		dNodes.resize( dNodes.size() + RES );
	}
	node->subnodes[ind] = &( dNodes[head] );
	omp_unset_lock( &lock );
	nNodes++;

	// Doesn't set occs as addRead...() will do that
	Node* tmpNode = node->subnodes[ind];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
	tmpNode->occs = 1;
	for( int i = 0; i < NBASES; i++ )
		tmpNode->subnodes[i] = nullptr;
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne( uint64_t readNum, short offset ) 
{
	//std::lock_guard<std::mutex> cncn(gtMut);
	std::vector<Node*> paths;

	Node* node = root;
	SeqRead* read = &GTH::seqReads[readNum];
	bool returnBool = 0, updateBool = 0;

	if( root->offset == offset && ( uint64_t )root->readNum == readNum )
		return;

	for( int i = offset + 1; i < GTH::seqReads[readNum].size(); i++ ) {
		short ind = ( *read ).getBaseInd( i );
		paths.push_back( node );

		omp_set_lock( &node->lock );
		if( !node->subnodes[ind] ) {
			//std::cout << ( *read ).getQual( i ) << std::endl;
			createNode( node, ind, ( *read ).getQual( i ), readNum, i );
			returnBool = updateBool = 1;
			if( countChildren( node ) == 1 )
				// Even if it doesn't balance it, a node was created and so
				// the path needs updating
				balanceNode( node );
		}
		if( i + 1 == GTH::seqReads[readNum].size() ) {
			// If the end is reached, they still count as occurrences of the path...
			updateBool = 1;
			// If a subnode exists and it wasn't created above
			if( node->subnodes[ind] && !returnBool)
				paths.push_back( node->subnodes[ind] );
			// For balancing purposes
			if( node->subnodes[ind] && !countChildren( node->subnodes[ind]) )
				potAddNode( node->subnodes[ind] );
		}
		omp_unset_lock( &node->lock );

		// It will always enter updateBool, just needs to know when
		if( updateBool )
			for( unsigned int j = 0, k = offset; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *read ).getQual(k) );
		if( returnBool )
			break;

		node = node->subnodes[ind];
	}
	paths.clear();
}

void GTree::potAddNode(Node* node)
{
	int64_t rN = node->readNum;
	SeqRead* tRead = &GTH::seqReads[rN];
	short offset = node->offset + 1;

	if( offset == ( *tRead ).size() )
		return;

	short ind = ( *tRead ).getBaseInd( offset );
	char qual = ( *tRead ).getQual( offset );

	createNode( node, ind, qual, rN, offset );
}

void GTree::balanceNode( Node* node )
{
	// Get the offset and read before overiding/updating
	int64_t lReadNum = node->readNum;
	SeqRead* lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	bool clearBool = 0;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if( lOffset == ( *lRead ).size() )
		return;

	short lInd = ( *lRead ).getBaseInd( lOffset );
	char lQual = ( *lRead ).getQual( lOffset );
	
	if(node->subnodes[lInd] && 
			lOffset == node->subnodes[lInd]->offset && 
			lReadNum == node->subnodes[lInd]->readNum)
		return;

	// If the paths are different:
	if( !node->subnodes[lInd] ) {
		createNode( node, lInd, lQual, lReadNum, lOffset );
		return;
	}
	node = node->subnodes[lInd];
	updateWeightAndOccs(node, lQual);

	// If the paths follow the same route:
	int64_t rReadNum = node->readNum;
	SeqRead* rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->offset + 1;
	lOffset++;

	std::vector<short> lIndPath;
	std::vector<short> rIndPath;
	short tmpLOffset = lOffset, tmpROffset = rOffset;

	while( lOffset < ( *lRead ).size() && rOffset < ( *rRead ).size() ) {
		lIndPath.push_back( ( *lRead ).getBaseInd( lOffset ) );
		rIndPath.push_back( ( *rRead ).getBaseInd( rOffset ) );
		if(( *rRead ).getBaseInd( rOffset ) != ( *lRead ).getBaseInd( lOffset )) {
			clearBool = 1;
			break;
		}
		lOffset++;
		rOffset++;
	} 

	// Update to the end of shared bases
	for( unsigned short i = 0; i < lIndPath.size() - clearBool; i++ ) {
		createNode( node, lIndPath[i], ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
		node = node->subnodes[lIndPath[i]];
		updateWeightAndOccs( node, ( *rRead ).getQual( tmpROffset ) );
		tmpLOffset++;
		tmpROffset++;
	}

	// If they differ, create a node each
	if( clearBool ) {
		createNode( node, lIndPath.back(), ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
		createNode( node, rIndPath.back(), ( *rRead ).getQual( tmpROffset ), rReadNum, tmpROffset );
		return;
	}

	// If one ends, create the relevant node
	if( lOffset < ( *lRead ).size() ) {
		lInd = ( *lRead ).getBaseInd( tmpLOffset );
		lQual = ( *lRead ).getQual( tmpLOffset );
		createNode( node, lInd, lQual, lReadNum, tmpLOffset );
	} 
	if( rOffset < ( *rRead ).size() ) {
		short rInd = ( *rRead ).getBaseInd( tmpROffset );
		char rQual = ( *rRead ).getQual(tmpROffset );
		createNode( node, rInd, rQual, rReadNum, tmpROffset );
	}

	return;
}
