#include "Node.H"

Node::Node(): occs(0), weight(0), offset(0), readNum(0)
{
	omp_init_lock(&lock);
	for(int i = 0; i < NBASES; i++)
		subnodes[i] = nullptr;
}

Node::Node(const Node& tmpNode) 
{
	float tmpWeight = tmpNode.weight;
	int64_t tmpOccs = tmpNode.occs;

	omp_init_lock(&lock);
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
	int64_t tmpOccs = tmpNode.occs;

	omp_init_lock(&lock);
	this->weight = tmpWeight;
	this->occs = tmpOccs;
	this->offset = tmpNode.offset;
	this->readNum = tmpNode.readNum;

	for(int i = 0; i < NBASES; i++)
		this->subnodes[i] = tmpNode.subnodes[i];

	return *this;
}

