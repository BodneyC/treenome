#include <sstream>
#include "includes/GTree.H"

namespace mingw_fix {
	template< typename T > std::string to_string(const T &val)
	{
		std::ostringstream ss;
		ss << val ;
		return ss.str() ;
	}
}

void GTree::addReadFULL(std::string &read)
{
	// Very primitive stage at the mo
	Node *tmpNode = root;
	root->occurences++;

	for(unsigned int i = 1; i < read.length(); i++) {
		short ind = (read[i] & 0xF) >> 1;
		if(!(tmpNode->subnodes[ind]))
			tmpNode->subnodes[ind] = new Node;
		tmpNode = tmpNode->subnodes[ind];
		tmpNode->occurences++;
	}
}

void GTree::addReadONE()
{

}

void GTree::deleteTree(Node* node)
{
	if(node) {
		for(int i = 0; i < 4; i++)
			deleteTree(node->subnodes[i]);
		delete node;
	}
}

// Definately needs work
void GTree::printAllPaths(Node *node, int len, short label)
{
    if(!node)
        return;

	occuPaths.erase(len, occuPaths.length());
	std::string val = mingw_fix::to_string(node->occurences);
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
	default:
		basePaths += 'N';
		break;
	}
	for(unsigned int i = 0; i < val.length(); i++)
		basePaths += '-';

	len += val.length() + 1;

	bool check = 0;
	for(int i = 0; i < 4; i++)
		if(node->subnodes[i])
			check = true;
			
    if(!check) {
		occuPaths.erase(occuPaths.length() - 1);
		std::cout << occuPaths << std::endl;
		basePaths.erase(basePaths.length() - 1);
		std::cout << basePaths << "\n" << std::endl;
        return;
    }

	for(int i = 0; i < 4; i++)
		printAllPaths(node->subnodes[i], len, i);
}
