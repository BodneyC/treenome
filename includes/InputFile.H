#include <fstream>
#include <vector>
#include <string>
#include "GTree.H"

#define FILE_ERROR -2

// Though sparse, this class may be expanded upon for processing/changing
// of storage mnethod
class InputFile {
public:
	long nReads, readLength;

	InputFile(const std::string &_filename):
		filename(_filename) {}

	bool readFastQ(
			std::vector<std::string> &reads,
			std::vector<std::string> &qualities
	);

private:
	std::string filename;
};
