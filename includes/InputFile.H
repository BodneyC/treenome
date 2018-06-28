#include <fstream>
#include <vector>
#include <string>
#include "GTree.H"

#define FILE_ERROR -1

class InputFile {
public:
	long nReads, readLength;

	InputFile(const std::string &_filename):
		filename(_filename) {}

	void readFastQ(TreeTop &treeTop);

private:
	std::string filename;
};
