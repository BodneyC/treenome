#include "includes/InputFile.H"

bool InputFile::readFastQ(
		std::vector<std::string> &reads,
		std::vector<std::string> &qualities
		)
{
	long i = 3, j = 1;
	std::ifstream inpFile(filename.c_str());
	if(!inpFile) {
		return 0;
	}

	for(std::string line; std::getline(inpFile, line); i++, j++) {
		if(!(i % 4))
			reads.push_back(line);
		if(!(j % 4))
			qualities.push_back(line);
	}

	nReads = (i - 3) / 4;
	readLength = reads[0].length();

	return 1;
}
