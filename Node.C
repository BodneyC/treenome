#include "includes/GTree.H"

Node::Node(): weight(0), offset(0), readNum(0) 
{
	occs = 0;
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = nullptr;
}

Node::Node(const Node& tmpNode) 
{
	float tmpWeight = tmpNode.weight;
	long tmpOccs = tmpNode.occs;

	weight = tmpWeight;
	occs = tmpOccs;
	offset = tmpNode.offset;
	readNum = tmpNode.readNum;
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = tmpNode.subnodes[i];
}

Node& Node::operator=(const Node& tmpNode) 
{
	float tmpWeight = tmpNode.weight;
	long tmpOccs = tmpNode.occs;

	this->weight = tmpWeight;
	this->occs = tmpOccs;
	this->offset = tmpNode.offset;
	this->readNum = tmpNode.readNum;

	for(int i = 0; i < NBASES; i++)
		this->subnodes[i] = tmpNode.subnodes[i];

	return *this;
}


LeafNode::LeafNode(const LeafNode& tmpLeafNode) 
{
	float tmpWeight = tmpLeafNode.weight;
	long tmpOccs = tmpLeafNode.occs;

	weight = tmpWeight;
	occs = tmpOccs;
	offset = tmpLeafNode.offset;
	readNum = tmpLeafNode.readNum;
}

LeafNode& LeafNode::operator=(const LeafNode& tmpLeafNode) 
{
	float tmpWeight = tmpLeafNode.weight;
	long tmpOccs = tmpLeafNode.occs;

	this->weight = tmpWeight;
	this->occs = tmpOccs;
	this->offset = tmpLeafNode.offset;
	this->readNum = tmpLeafNode.readNum;

	return *this;
}

