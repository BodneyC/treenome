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

char retLabel(int label)
{
	switch(label) {
	case 0:
		return 'A'; break;
	case 1:
		return 'C'; break;
	case 2:
		return 'T'; break;
	case 3:
		return 'G'; break;
	default:
		return 'N'; break;
	}
}

float cumAve(float curAve, char qual, ulong occ)
{
	// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
	return curAve + 
		((static_cast<float>(qual) - curAve) / static_cast<double>(occ));
}

/** ---------------- Graph Creation ---------------- **/
void GTree::createRoot(std::vector<std::string> &reads, int label)
{
	root = new Node;
	char lab = retLabel(label);

	for(uint i = 0; i < reads.size(); i++) {
		auto offset = reads[i].find(lab, 0);
		if(offset != std::string::npos) {
			root->readNum = i;
			root->offset = offset;
			return;
		}
	}
}

void GTree::addReadOne(std::string &read, std::string &qual)
{

}

void GTree::addReadFull(ulong readNum, short offset, std::string &read, std::string &qual)
{
	Node *tmpNode = root;
	tmpNode->occurences++;
	tmpNode->weight = cumAve(root->weight, qual[0], tmpNode->occurences);
	
	for(uint i = offset + 1; i < read.length(); i++) {
		short ind = (read[i] & 0xF) >> 1;
		if(!(tmpNode->subnodes[ind])) {
			tmpNode->subnodes[ind] = new Node;
			tmpNode->subnodes[ind]->readNum = readNum;
			tmpNode->subnodes[ind]->offset = offset;
		}
		tmpNode = tmpNode->subnodes[ind];
		tmpNode->occurences++;
		tmpNode->weight = cumAve(tmpNode->weight, qual[0], tmpNode->occurences);
	}
}

/** ------------- None recursive calls ------------- **/
void GTree::cleanBranchesNR(std::string &read)
{
	Node *tmpNode = root;

	for(uint i = 1; i < read.length(); i++) {
		short ind = (read[i] & 0xF) >> 1;
		if(tmpNode->occurences == 1 && 
				tmpNode->subnodes[ind]->occurences == 1) {
			deleteTreeNR(&tmpNode->subnodes[ind], read.substr(i, read.length()));
			return;
		}
		tmpNode = tmpNode->subnodes[ind];
	}
}

// If both occurences are both 1, the remaining branch is a singly linked line
void GTree::deleteTreeNR(Node **node, std::string readSub)
{
	Node **curNode = node;
	
	for(uint i = 0; i < readSub.length(); i++) {
		short ind = (readSub[i + 1] & 0xF) >> 1;

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
	if(nodeCount == 1 && node->occurences == 1 && 
			node->subnodes[ind]->occurences == 1) {
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
	std::string val = mingw_fix::to_string(node->occurences);
#else
	std::string val = std::to_string(node->occurences);
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
