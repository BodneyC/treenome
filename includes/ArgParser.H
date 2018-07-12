/********************************************************************
 * Filename: ArgParser.H [C++ header code]
 *
 * Description: Simple class for command line argument parsing
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _ARGPARSER_
#define _ARGPARSER_

#include <string>
#include <vector>
#include <fstream>

#define USAGE_ERROR -1
#define FILE_ERROR -2

struct CMDArgs {
	std::string iFilename;
	bool printToScreen;
	bool preProcess;

	CMDArgs(): iFilename(""), printToScreen(0), preProcess(0) {  }
};

class ArgParser {
public:
	ArgParser(int& argc, char** argv) 
	{
		this->argc = argc;
		for(int i = 1; i < argc; i++)
			args.push_back(std::string(argv[i]));
		this->argc--;
	}
	int fillArgs(struct CMDArgs &argList)
	{
		for(int i = 0; i < argc; i++) {
			if(args[i] == "-f") {
				i++;
				std::ifstream testFile(args[i]);
				if(testFile.good())
					argList.iFilename = args[i];
				else 
					return FILE_ERROR;
			} else if (args[i] == "-o"){
				argList.printToScreen = 1;
			} else if(args[i] == "-p") {
				argList.preProcess = 1;
			} else {
				return USAGE_ERROR;
			}
		}
		return 0;
	}
	bool argExists(const std::string &arg)
	{
		for(int i = 0; i < argc; i++)
			if(args[i] == arg)
				return true;
		return false;
	}
	const std::string stringOption(const std::string &arg)
	{
		for(int i = 0; i < argc; i++) 
			if(args[i] == arg && i + 1 < argc)
				return args[i + 1];
		return "NONE";
	}
private:
	int argc;
	std::vector<std::string> args;
};

#endif /*_ARGPARSER_*/
