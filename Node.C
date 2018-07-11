#include "includes/GTree.H"

Node::Node(): occs(0), weight(0), offset(0), readNum(0) 
{
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = nullptr;
}

Node::Node(const Node& tmpNode) 
{
	occs = weight = offset = readNum = 0;
	std::atomic<Node*> tmpAtomNode;
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = tmpAtomNode.load();
}

Node& Node::operator=(const Node& tmpNode) 
{
	this->occs = tmpNode.occs;
	this->weight = tmpNode.weight;
	this->offset = tmpNode.offset;
	this->readNum = tmpNode.readNum;

	for(int i = 0; i < NBASES; i++)
		this->subnodes[i] = tmpNode.subnodes[i].load();

	return *this;
}

