#pragma once

#include <map>

class Parser {

	const int argc;
	char * const *argv;
	std::map<std::string, std::string> parameters;

	void add_parameter_from_stream(const char *key, const std::stringstream& value);

public:
	Parser(const int argc, char * const *argv);
	void parse();
	friend std::ostream& operator<<(std::ostream& os, const Parser& parser);

};
