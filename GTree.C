#include "includes/GTree.H"

/** --------------- Helper Functions --------------- **/
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

namespace GTreeFuncs {
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
			return 'N';
		}
	}

	short mostOccs(Node *node)
	{
		short ind;
		ulong occCnt = 0;
		for(int i = 0; i < NBASES; i++) {
			if(node->subnodes[i] && node->subnodes[i]->occs > occCnt) {
				occCnt = node->subnodes[i]->occs;
				ind = i;
			}
		}
		return ind;
	}
}

void GTree::updateWeight(Node *node, char qual)
{
	// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
	node->weight += ((static_cast<float>(qual) - node->weight) / 
			static_cast<double>(++node->occs));
	//return curAve + 
	//	((static_cast<float>(qual) - curAve) / static_cast<double>(occ));
}

short GTree::countChildren(Node *node)
{
	short children = 0;

	for(int i = 0; i < NBASES; i++)
		if(node->subnodes[i])
			children++;

	return children;
}

/** --------------- Genome Creation ---------------- **/
void GTree::followPath(Node *node, short ind, std::string &sequence)
{
	short children;
	do {
		children = countChildren(node);
		sequence += GTreeFuncs::retLabel(ind);
		ind = GTreeFuncs::mostOccs(node);
		node = node->subnodes[ind];
	} while(children);
}

void GTree::addToSeq(ulong offset, std::string &sequence)
{
	Node *node = root;
	short ind = 0;
	uint i, j, seqLength = sequence.length() - offset;
	Node **path = new Node*[seqLength];

	for(i = 0; i < sequence.length() - offset; i++)
		path[i] = nullptr;

	// End of current sequence	
	for(i = offset + 1, j = 0; i < sequence.length(); i++, j++) {
		ind = BASE_IND(sequence[i]);
		if(node->subnodes[ind] && countChildren(node->subnodes[ind])) {
			path[j] = node;
			node = node->subnodes[ind];
		} else {
			delete[] path;
			return;
		}
	}

	ind = BASE_IND(sequence[i]);
	while(countChildren(node)) {
		ind = GTreeFuncs::mostOccs(node);
		sequence += GTreeFuncs::retLabel(ind);
		node->occs--;
		if(!node->subnodes[ind])
			break;
		node = node->subnodes[ind];
	}
	
	for(i = 0; i < seqLength; i++)
		if(path[i]) 
			path[i]->occs--;

	delete[] path;
}

/** ---------------- Graph Creation ---------------- **/
void GTree::createRoot(short ind, VecString &reads, VecString &quals)
{
	root = new Node;
	char lab = GTreeFuncs::retLabel(ind);

	for(uint i = 0; i < reads.size(); i++) {
		auto offset = reads[i].find(lab, 0);
		if(offset != std::string::npos) {
			char qual = quals[i][offset];
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

	// addReadFull() looses accuracy with this
	// updateWeight(node, qual);
}

void GTree::balanceNode(Node *node, VecString &reads, VecString &quals)
{
	// Get the offset and read before overiding/updating
	ulong lReadNum = node->readNum;
	std::string lRead = reads[lReadNum];
	short lOffset = node->offset + 1;
	short lInd = BASE_IND(lRead[lOffset]);
	
	// If the paths are different:
	if(!node->subnodes[lInd]) {
		createNode(node, lInd, quals[lReadNum][lOffset]);
		updateWeight(node->subnodes[lInd], quals[lReadNum][lOffset]);
		return;
	}

	// If the paths are literally the same:
	if(lOffset == node->subnodes[lInd]->offset &&
			lReadNum == node->subnodes[lInd]->readNum)
		return;

	// If the paths follow the same route:
	ulong rReadNum = node->subnodes[lInd]->readNum;
	std::string rRead = reads[rReadNum];
	short rOffset = rRead[node->subnodes[lInd]->offset + 1];
	short rInd = 0;
	lOffset++;
	while(lOffset < (short)lRead.length() && rOffset < (short)rRead.length()) {
		node = node->subnodes[lInd];
		updateWeight(node, quals[rReadNum][rOffset]);
		rInd = BASE_IND(rRead[rOffset]);
		lInd = BASE_IND(lRead[lOffset]);

		if(lRead[lOffset] != rRead[rOffset]) {
			createNode(node, lInd, quals[lReadNum][lOffset]);
			updateWeight(node->subnodes[lInd], quals[lReadNum][lOffset]);
			createNode(node, rInd, quals[rReadNum][rOffset]);
			updateWeight(node->subnodes[rInd], quals[rReadNum][rOffset]);
			return;
		}

		lOffset++;
		rOffset++;
	} 

	// There will be imbalances caused by running out of read length
	if(lOffset < (short)lRead.length()) {
		createNode(node, lInd, quals[lReadNum][lOffset]);
		updateWeight(node->subnodes[lInd], quals[lReadNum][lOffset]);
	} else {
		createNode(node, rInd, quals[rReadNum][rOffset]);
		updateWeight(node->subnodes[rInd], quals[rReadNum][rOffset]);
	}
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne(ulong readNum, short offset, 
		VecString &reads, VecString &quals)
{
	Node *node = root;
	
	for(uint i = offset + 1; i < reads[readNum].length(); i++) {
		short ind = BASE_IND(reads[readNum][i]);
		updateWeight(node, quals[readNum][i - 1]);
		if(!(node->subnodes[ind])) {
			createNode(node, ind, quals[readNum][i]);
			updateWeight(node->subnodes[ind], quals[readNum][i]);
			if(countChildren(node) == 1)
				balanceNode(node, reads, quals);
			return;
		}
		node = node->subnodes[ind];
	}
}

void GTree::addReadFull(ulong readNum, short offset, 
		std::string &read, 
		std::string &qual)
{
	Node *node = root;
	
	for(uint i = offset + 1; i < read.length(); i++) {
		updateWeight(node, qual[i - 1]);
		short ind = BASE_IND(read[i]);
		if(!(node->subnodes[ind]))
			createNode(node, ind, qual[i]);
		node = node->subnodes[ind];
	}
}

/** ------------- None recursive calls ------------- **/
Node* GTree::cleanBranchesNR(short offset, std::string &read)
{
	Node *nextNode = root, *curNode;

	// Follow path in search of branch of ones (1->1->1->1...)
	for(uint i = offset + 1; i < read.length(); i++) {
		short ind = BASE_IND(read[i]);
		if(nextNode->occs == 1 && 
				nextNode->subnodes[ind]->occs == 1) {
			// Delete branch of ones leaving single unique path
			deleteTreeLL(&nextNode->subnodes[ind], 
					read.substr(i, read.length()));
			// Count number of children of parent, 
			//	if only 1 then node need balancing
			if(countChildren(curNode) < 2) {
				return curNode;
			}
			break;
		}
		// Traverse path of read maintaining previous node of occs > 1
		curNode = nextNode;
		nextNode = nextNode->subnodes[ind];
	}
	return nullptr;
}

// If both occs are both 1, the remaining branch is a simple linked line
void GTree::deleteTreeLL(Node **node, std::string readSub)
{
	Node **curNode = node;
	
	for(uint i = 0; i < readSub.length(); i++) {
		short ind = BASE_IND(readSub[i + 1]);
		Node *nextNode = (*curNode)->subnodes[ind];

		delete *(curNode);
		*curNode = nullptr;
		*curNode = nextNode;
	}
	*curNode = nullptr;
}

/** -------------- Recursive Cleaning -------------- **/
void GTree::cleanBranches(Node *node)
{
	if(!node)
		return;

	int nodeCount = 0;
	int ind = 0;
	for(int i = 0; i < NBASES; i++) {
		if(node->subnodes[i]) {
			nodeCount++;
			ind = i;
		}
		if(nodeCount == 2)
			break;
	}
	if(nodeCount == 1 && node->occs == 1 && 
			node->subnodes[ind]->occs == 1) {
		deleteTree(node->subnodes[ind]);
		node->subnodes[ind] = nullptr;
		return;
	} 
	for(int i = 0; i < NBASES; i++)
		cleanBranches(node->subnodes[i]);
}

void GTree::deleteTree(Node* node)
{
	if(node) {
		for(int i = 0; i < NBASES; i++)
			deleteTree(node->subnodes[i]);
		delete node;
		node = nullptr;
	}
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
	switch (label) {
	case 0:
		basePaths += 'A';
		break;
	case 1:
		basePaths += 'C';
		break;
	case 2:
		basePaths += 'T';
		break;
	case 3:
		basePaths += 'G';
		break;
	default: // Shouldn't occur
		basePaths += 'N';
		break;
	}
	for(uint i = 0; i < val.length(); i++)
		basePaths += '-';
	len += val.length() + 1;
	bool check = 0;
	for(int i = 0; i < NBASES; i++)
		if(node->subnodes[i])
			check = true;
    if(!check) {
		occuPaths.erase(occuPaths.length() - 1);
		std::cout << occuPaths << std::endl;
		basePaths.erase(basePaths.length() - 1);
		std::cout << basePaths << "\n" << std::endl;
        return;
    }
	for(int i = 0; i < NBASES; i++)
		printAllPaths(node->subnodes[i], len, i);
}
