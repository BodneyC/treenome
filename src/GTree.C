/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree class, templatized for 
 *		Node/pNode purposes
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/GTree.H"

#ifdef __MINGW64__
namespace mingw_fix {
	{
		std::ostringstream ss;
		ss << val ;
		return ss.str() ;
	}
}
#endif /*__MINGW32__*/

/** --------------- Helper Functions --------------- **/
namespace GTH {
	std::vector<SeqRead> seqReads;
	double thresh;
	double startWeights[NBASES];
	int64_t startOccs[NBASES];
	int scoreSys;
	int32_t readLength;

	char retLabel( int label )
	{
		switch( label ) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'T';
		case 3:
			return 'G';
		default: // 7 but it shouldn't
			return 'P';
		}
	}

	void removeDoubleEnding( std::string& doubleString )
	{
		doubleString.erase( doubleString.find_last_not_of('0') + 1, std::string::npos );
		doubleString.erase( doubleString.find_last_not_of('.') + 1, std::string::npos );
	}
	
	template <typename T>
	std::string valToString( T& val )
	{
#ifdef __MINGW64__
		return mingw_fix::to_string( val );
#else
		return std::to_string( val );
#endif /*__MINGW64__*/
	}
}

/** ------------- - Helper Functions --------------- **/
Node* GTree::getRoot()
{
	if( root ) 
		return root; 
	else 
		return nullptr; 
}

int64_t GTree::getNNodes()
{
	return nNodes;
}

short GTree::countChildren( Node* node )
{
	short children = 0;

	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			children++;

	return children;
}

signed short GTree::highestThresh( Node* node )
{
	short ind = -1;
	double maxRat = GTH::thresh;

	for( short i = 0; i < NBASES; i++ ) {
		Node* tmpNode = &dNodes[node->subnodes[i]];
		if( node->subnodes[i] && tmpNode->occs ) {
			double curRat = tmpNode->getRatio();
			if( curRat > maxRat ) {
				maxRat = curRat;
				ind = i;
			}
		}
	}

	return ind;
}

/** --------------- Genome Creation ---------------- **/
void GTree::reduceOccsAndWeight( std::vector<int32_t>& paths )
{
	for( uint32_t i = 0; i < paths.size(); i++ ) {
		dNodes[paths[i]].occs--;
		dNodes[paths[i]].weight = dNodes[paths[i]].weight.load() - 1;
	}
}

Node* GTree::followSeq( int32_t offset, std::string& sequence, std::vector<int32_t>& paths )
{
	Node* node = root;
	paths.push_back( 0 );

	// Reach the end of the existing path
	for( uint32_t i = offset; i < sequence.length(); i++ ) {
		short ind = BASE_IND( sequence[i] );
		if( node->subnodes[ind] ) {
			paths.push_back( node->subnodes[ind] );
			node = &dNodes[node->subnodes[ind]];
		} else {
			return nullptr;
		}
	}

	if( highestThresh( node ) == -1 )
		return nullptr;

	return node;
}

void GTree::followBranch( Node* node, short ind, std::string &sequence )
{
	short children;
	do {
		children = countChildren( node );
		sequence += GTH::retLabel( ind );
		ind = highestThresh( node );
		node->occs--;
		node->weight = node->weight.load() - 1;
		if( ind == -1 )
			return;
		node = &dNodes[node->subnodes[ind]];
	} while( children );
}

void GTree::addToSeq( uint64_t offset, std::string &sequence )
{
	Node* node = root;
	short ind = 0;
	std::vector< Node* > nPath;

	// End of current sequence	
	for( uint32_t i = offset + 1; i < sequence.length(); i++ ) {
		ind = BASE_IND( sequence[i] );
		if( node->subnodes[ind] && countChildren( &dNodes[node->subnodes[ind]] ) ) {
			nPath.push_back( node );
			node = &dNodes[node->subnodes[ind]];
		} else {
			return;
		}
	}

	ind = highestThresh( node );
	if( ind != -1 ) {
		followBranch( &dNodes[node->subnodes[ind]], ind, sequence );
		
		// Only if something is contributed to the sequence should the 
		// occurences be lowered
		for( uint32_t i = 0; i < nPath.size(); i++ ) {
			nPath[i]->occs--;
			nPath[i]->weight = nPath[i]->weight.load() - 1;
		}
	}
}

/** --------------- Tree from Reads ---------------- **/
void GTree::init( short i ) {
	// A deque has been used to prevent frequent reallocations of existing
	// vector data when dNodes is resized
	omp_init_lock( &lock );
	dNodes.resize( RES );
	createRoot( i );
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
	char qual;
	root = &( dNodes[0] );
	nNodes++;

	for( int64_t i = 0; i < GTH::seqReads.size(); i++ ) {
		int offset = -1;
		// Find algorithm in bool vec
		for( int64_t j = 0; j < GTH::seqReads[i].size(); j++ ) {
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
	if( ( unsigned ) head == dNodes.size() ) {
		dNodes.resize( dNodes.size() + RES );
	}
	node->subnodes[ind] = head;
	omp_unset_lock( &lock );
	nNodes++;

	// Doesn't set occs as addRead...() will do that
	Node* tmpNode = &dNodes[node->subnodes[ind]];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
	tmpNode->occs = 1;
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne( uint64_t readNum, short offset ) 
{
	std::vector<Node*> paths;

	Node* node = root;
	SeqRead* read = &GTH::seqReads[readNum];
	bool returnBool = 0, updateBool = 0;

	// Edge case
	if( root->readNum == readNum && root->offset == offset )
		return;

	for( int i = offset + 1; i < GTH::seqReads[readNum].size(); i++ ) {
		short ind = ( *read ).getBaseInd( i );
		paths.push_back( node );

		omp_set_lock( &node->lock );
		int64_t nextNodeLoc = node->subnodes[ind];
		if( !nextNodeLoc ) {
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
			if( nextNodeLoc && !returnBool)
				paths.push_back( &dNodes[nextNodeLoc] );
			// For balancing purposes
			if( nextNodeLoc && !countChildren( &dNodes[nextNodeLoc]) )
				potAddNode( &dNodes[nextNodeLoc] );
		}
		omp_unset_lock( &node->lock );

		// It will always enter updateBool, just needs to know when
		if( updateBool )
			for( unsigned int j = 0, k = offset; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *read ).getQual(k) );
		if( returnBool )
			break;

		node = &dNodes[nextNodeLoc];
	}
	paths.clear();
}

void GTree::balanceNode( Node* node )
{
	// Get the offset and read before overiding/updating
	uint64_t lReadNum = node->readNum;
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
			lOffset == dNodes[node->subnodes[lInd]].offset && 
			lReadNum == dNodes[node->subnodes[lInd]].readNum)
		return;

	// If the paths are different:
	if( !node->subnodes[lInd] ) {
		createNode( node, lInd, lQual, lReadNum, lOffset );
		return;
	}
	node = &dNodes[node->subnodes[lInd]];
	updateWeightAndOccs(node, lQual);

	// If the paths follow the same route:
	uint64_t rReadNum = node->readNum;
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
		node = &dNodes[node->subnodes[lIndPath[i]]];
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

/** ---------------- Tree Storage ------------------ **/
void GTree::writeTreeToFile( std::ofstream& storeFile )
{
	for( uint64_t i = 0; i < dNodes.size(); i++ )
		storeFile.write( ( char* )&dNodes[i], sizeof( Node ) );
}

/** --------------- Tree From File ----------------- **/
void GTree::resizeVector( int64_t tmp64 )
{
	nNodes = tmp64;
	if( tmp64 % RES )
		tmp64 = (tmp64 / RES + 1) * RES;
	//std::cout << tmp64 << std::endl;
	dNodes.resize( tmp64 );
}

void GTree::writeVector( std::ifstream& inFile )
{
	for( int64_t i = 0; i < dNodes.size(); i++ ) {
		Node tmpNode;
		inFile.read( ( char* ) &tmpNode, sizeof( Node ) );
		dNodes[i] = tmpNode;
	}
	root = &dNodes[0];
}

/** -------------- Useful Functions ---------------- **/
int32_t GTree::maxDepth( Node* node )
{
	int32_t depth[NBASES] = { -1 };
	int32_t ret = 0;

	for( int i = 0; i < NBASES; i++ ) {
		if( node->subnodes[i] ) {
			ret = -1;
			depth[i] = maxDepth( &dNodes[node->subnodes[i]] );
		}
	}

	if( !ret )
		return 1;
	
	for( int i = 0; i < NBASES; i++ )
		if( depth[i] > ret )
			ret = depth[i];

	return ret + 1;
}

int32_t GTree::maxDepthAThresh( Node* node )
{
	int32_t depth[NBASES] = { -1 };
	int32_t ret = 0;

	if( node->getRatio() < GTH::thresh )
		return 0;

	for( int i = 0; i < NBASES; i++ ) {
		if( node->subnodes[i] ) {
			ret = -1;
			depth[i] = maxDepthAThresh( &dNodes[node->subnodes[i]] );
		}
	}

	if( !ret )
		return 1;
	
	for( int i = 0; i < NBASES; i++ )
		if( depth[i] > ret )
			ret = depth[i];

	return ret + 1;
}

void GTree::printAllPaths( short label ) { 
	printAllPaths( root, 0, label );
	basePaths = occuPaths = "";
}

void GTree::printAllPaths( Node* node, int len, short label )
{
	occuPaths.erase( len, occuPaths.length() );
	std::string val = GTH::valToString( node->occs );
    occuPaths += val;
	occuPaths += "-";
	basePaths.erase( len, basePaths.length() );
	basePaths += GTH::retLabel( label );
	for( uint32_t i = 0; i < val.length(); i++ )
		basePaths += '-';
	len += val.length() + 1;
	bool check = 0;
	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			check = true;
	//std::cout << "WEIGHT: " << node->weight << '\n';
    if( !check ) {
		occuPaths.erase( occuPaths.length() - 1 );
		std::cout << occuPaths << ": EOS" << std::endl;
		basePaths.erase( basePaths.length() - 1 );
		std::cout << basePaths << ": EOS\n" << std::endl;
        return;
    }
	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			printAllPaths( &dNodes[node->subnodes[i]], len, i );
}

