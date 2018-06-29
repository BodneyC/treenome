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

	bool readFastQ(
			std::vector<std::string> &reads,
			std::vector<std::string> &qualities
			);

private:
	std::string filename;
};
