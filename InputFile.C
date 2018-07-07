#include "includes/InputFile.H"

bool InputFile::readFastQ()
{
	std::string seqLine = "";
	std::string qualLine = "";
	long i = 3, j = 1;
	short len = 0;

	std::ifstream inpFile(filename.c_str());
	if(!inpFile) {
		return 0;
	}

	// Bit of a hacky way of doing this...
	for(std::string line; std::getline(inpFile, line); i++, j++) {
		if(i == 3)
			readLength = line.length();
		if(!(i % 4)) {
			len = line.find('N');
			line = line.substr(0, len);

			seqLine = line;
		}
		if(!(j % 4)) {
			line = line.substr(0, len);

			qualLine = line;
			GTH::seqReads.push_back(SeqRead(seqLine, qualLine));
		}
	}

	nReads = (i - 3) / 4;

	return 1;
}
