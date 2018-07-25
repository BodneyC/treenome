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
#define IN_FILE_ERROR -2
#define OUT_FILE_ERROR -3

struct CMDArgs {
	std::string iFilename, oFilename, lFilename;
	bool printToScreen, storeToFile, loadFile;

	CMDArgs(): 
		iFilename(""), oFilename(""), lFilename(""),
		printToScreen(0), storeToFile(0), loadFile(0) {  }
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
				if(inFileCheck(args[i]))
					return IN_FILE_ERROR;
				argList.iFilename = args[i];
			} else if(args[i] == "-s") {
				i++;
				if(outFileCheck(args[i]))
					return OUT_FILE_ERROR;
				argList.oFilename = args[i];
				argList.storeToFile = 1;
			} else if (args[i] == "-l") {
				i++;
				if(inFileCheck(args[i]))
					return IN_FILE_ERROR;
				argList.lFilename = args[i];
				argList.loadFile = 1;
			} else if(args[i] == "-o"){
				argList.printToScreen = 1;
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
	int inFileCheck(std::string filename)
	{
		std::ifstream testFile(filename);
		if(testFile.good())
			return 0;
		else
			return IN_FILE_ERROR;
	}
	int outFileCheck(std::string filename)
	{
		std::ofstream testFile(filename);
		if(testFile.good())
			return 0;
		else
			return OUT_FILE_ERROR;
	}
private:
	int argc;
	std::vector<std::string> args;
};

#endif /*_ARGPARSER_*/
