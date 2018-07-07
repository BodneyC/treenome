#include "includes/InputFile.H"

bool InputFile::readFastQ(
		std::vector<std::string> &reads,
		std::vector<std::string> &qualities)
{
	std::string seqLine = "";
	std::string qualLine = "";
	long i = 3, j = 1;
	short len = 0;

	std::ifstream inpFile(filename.c_str());
	if(!inpFile) {
		return 0;
	}

	for(std::string line; std::getline(inpFile, line); i++, j++) {
		if(i == 3)
			readLength = line.length();
		if(!(i % 4)) {
			len = line.find('N');
			line = line.substr(0, len);
			reads.push_back(line);

			seqLine = line;
		}
		if(!(j % 4)) {
			line = line.substr(0, len);
			qualities.push_back(line);

			qualLine = line;
			GTH::seqReads.push_back(SeqRead(seqLine, qualLine));
		}
	}

	nReads = (i - 3) / 4;

	return 1;
}
