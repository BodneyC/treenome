#ifndef _ARGPARSER_
#define _ARGPARSER_

#include <string>
#include <vector>

class ArgParser {
public:
	ArgParser(int& argc, char** argv) 
	{
		this->argc = argc;
		for(int i = 0; i < argc; i++)
			args.push_back(std::string(argv[i]));
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
