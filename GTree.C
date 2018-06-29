#include <sstream>
#include "includes/GTree.H"

#ifdef __MINGW64__
namespace mingw_fix {
	template< typename T > std::string to_string(const T &val)
	{
		std::ostringstream ss;
		ss << val ;
		return ss.str() ;
	}
}
#endif /*__MINGW32__*/

#ifdef READSONE
void GTree::addReadOne()
{

}
#endif /*READSONE*/

void GTree::addReadFull(std::string &read)
{
	// Very primitive stage at the mo
	Node *tmpNode = root;
	root->occurences++;

	for(uint i = 1; i < read.length(); i++) {
		short ind = (read[i] & 0xF) >> 1;
		if(!(tmpNode->subnodes[ind]))
			tmpNode->subnodes[ind] = new Node;
		tmpNode = tmpNode->subnodes[ind];
		tmpNode->occurences++;
	}
}

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

// Definately needs work
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
