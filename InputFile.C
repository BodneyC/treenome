#include "includes/InputFile.H"

void InputFile::readFastQ(TreeTop &treeTop)
{
	long i = 3, j = 1;
	std::ifstream inpFile(filename.c_str());
	if(!inpFile) {
		nReads = FILE_ERROR;
		return;
	}

	for(std::string line; std::getline(inpFile, line); i++, j++) {
		if(!(i % 4))
			treeTop.reads.push_back(line);
		if(!(j % 4))
			treeTop.qualities.push_back(line);
	}

	treeTop.nReads = nReads = (i - 3) / 4;
	treeTop.readLength = readLength = treeTop.reads[0].length();
}
