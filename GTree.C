/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "includes/GTree.H"

#ifdef __MINGW64__
#include <sstream>
namespace mingw_fix {
	template< typename T > std::string to_string(const T &val)
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

	char retLabel(int label)
	{
		switch(label) {
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

	short mostOccs(Node *node)
	{
		short ind = -1;
		long occCnt = 0;
		for(int i = 0; i < NBASES; i++) {
			if(node->subnodes[i] && node->subnodes[i]->occs > occCnt) {
				occCnt = node->subnodes[i]->occs;
				ind = i;
			}
		}
		return ind;
	}

	short countChildren(Node *node)
	{
		short children = 0;

		for(int i = 0; i < NBASES; i++)
			if(node->subnodes[i])
				children++;

		return children;
	}

	void updateWeight(Node *node, char qual)
	{
		// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
		node->weight += ((static_cast<float>(qual) - node->weight) / 
				static_cast<double>(++node->occs));
	}
}

/** --------------- Genome Creation ---------------- **/
void GTree::followPath(Node *node, short ind, std::string &sequence)
{
	short children;
	do {
		children = GTH::countChildren(node);
		sequence += GTH::retLabel(ind);
		ind = GTH::mostOccs(node);
		node->occs--;
		if(ind == -1)
			return;
		node = node->subnodes[ind];
	} while(children);
}

void GTree::addToSeq(long offset, std::string &sequence)
{
	Node *node = root;
	short ind = 0;
	uint i, j, seqLength = sequence.length() - offset;
	Node **path = new Node*[seqLength];

	for(i = 0; i < seqLength; i++)
		path[i] = nullptr;

	// End of current sequence	
	for(i = offset + 1, j = 0; i < sequence.length(); i++, j++) {
		ind = BASE_IND(sequence[i]);
		if(node->subnodes[ind] && GTH::countChildren(node->subnodes[ind])) {
			path[j] = node;
			node = node->subnodes[ind];
		} else {
			delete[] path;
			return;
		}
	}

	ind = GTH::mostOccs(node);
	if(ind != -1) {
		followPath(node->subnodes[ind], ind, sequence);
		
		// Only if something is contributed to the sequence should the 
		// occurences be lowered
		for(i = 0; i < seqLength; i++)
			if(path[i]) 
				path[i]->occs--;
	}

	delete[] path;
}

/** ---------------- Tree Creation ----------------- **/
void GTree::deleteTree(Node* node)
{
	if(node) {
		for(int i = 0; i < NBASES; i++)
			deleteTree(node->subnodes[i]);
		delete node;
		node = nullptr;
	}
}

void GTree::createRoot(short ind)
{
	root = new Node;

	for(uint i = 0; i < GTH::seqReads.size(); i++) {
		int offset = -1;
		// Find algorithm
		for(int j = 0; j < GTH::seqReads[i].size(); j++) {
			if(GTH::seqReads[i].getBaseInd(j) == ind) {
				offset = j;
				break;
			}
		}
		if(offset != -1) {
			char qual = GTH::seqReads[i].getQual(offset);
			root->readNum = i;
			root->offset = offset;
			root->weight = qual;
			// Doesn't set occs as addRead...() will do that
			return;
		}
	}
}

void GTree::createNode(Node *node, short ind, char qual)
{
	node->subnodes[ind] = new Node;
	Node *tmpNode = node->subnodes[ind];

	// Doesn't set occs as addRead...() will do that
	tmpNode->readNum = node->readNum;
	tmpNode->offset = node->offset + 1;
	tmpNode->weight = qual;

	GTH::updateWeight(tmpNode, qual);
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne(long readNum, short offset) 
{
	Node *node = root;
	SeqRead *read = &GTH::seqReads[readNum];
	
	for(int i = offset + 1; i < GTH::seqReads[readNum].size(); i++) {
		short ind = (*read).getBaseInd(i);
		GTH::updateWeight(node, (*read).getQual(i - 1));
		if(!(node->subnodes[ind])) {
			createNode(node, ind, (*read).getQual(i));
			if(GTH::countChildren(node) == 1)
				balanceNode(node);
			return;
		}
		node = node->subnodes[ind];
	}
}

void GTree::balanceNode(Node *node)
{
	// Get the offset and read before overiding/updating
	long lReadNum = node->readNum;
	SeqRead *lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	short lInd = (*lRead).getBaseInd(lOffset);
	char lQual = (*lRead).getQual(lOffset);
	
	// If the paths are different:
	if(!node->subnodes[lInd]) {
		createNode(node, lInd, lQual);
		return;
	}

	// If the paths are literally the same:
	if(lOffset == node->subnodes[lInd]->offset &&
			lReadNum == node->subnodes[lInd]->readNum)
		return;

	// If the paths follow the same route:
	long rReadNum = node->subnodes[lInd]->readNum;
	SeqRead *rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->subnodes[lInd]->offset + 1;
	short rInd = 0;
	char rQual;
	lOffset++;
	while(lOffset < (*lRead).size() && rOffset < (*rRead).size()) {
		lInd = (*lRead).getBaseInd(lOffset);
		rInd = (*rRead).getBaseInd(rOffset);
		lQual = (*lRead).getQual(lOffset);
		rQual = (*rRead).getQual(rOffset);

		node = node->subnodes[lInd];
		GTH::updateWeight(node, lQual);
		GTH::updateWeight(node, rQual);

		if((*rRead).getBaseInd(rOffset) != (*lRead).getBaseInd(lOffset)) {
			createNode(node, lInd, lQual);
			createNode(node, rInd, rQual);
			return;
		}

		lOffset++;
		rOffset++;
	} 

	// There will be imbalances caused by running out of read length
	lQual = (*lRead).getQual(lOffset);
	rQual = (*rRead).getQual(rOffset);
	if(lOffset < (*lRead).size())
		createNode(node, lInd, lQual);
	else
		createNode(node, rInd, rQual);
}

/** ---------------- Path Printing ----------------- **/
void GTree::printAllPaths(Node *node, int len, short label)
{
    if(!node)
        return;

	occuPaths.erase(len, occuPaths.length());
#ifdef __MINGW64__
	std::string val = mingw_fix::to_string(node->occs);
#else
	std::string val = std::to_string(node->occs);
#endif /*__MINGW32__*/
    occuPaths += val;
	occuPaths += "-";
	basePaths.erase(len, basePaths.length());
	basePaths += GTH::retLabel(label);
	for(uint i = 0; i < val.length(); i++)
		basePaths += '-';
	len += val.length() + 1;
	bool check = 0;
	for(int i = 0; i < NBASES; i++)
		if(node->subnodes[i])
			check = true;
    if(!check) {
		occuPaths.erase(occuPaths.length() - 1);
		std::cout << occuPaths << ": EOS" << std::endl;
		basePaths.erase(basePaths.length() - 1);
		std::cout << basePaths << ": EOS\n" << std::endl;
        return;
    }
	for(int i = 0; i < NBASES; i++)
		printAllPaths(node->subnodes[i], len, i);
}
