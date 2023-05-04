#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include "Parser.hpp"

Parser::Parser(const int argc, char * const *argv) : argc{argc}, argv{argv} {}

void Parser::parse()
{
	if (argc == 1)
		return;

	int i = 1;
	char *cur_key = nullptr;
	std::stringstream value;
	while (i < argc) {
		if (argv[i][0] == '-') {
			if (cur_key != nullptr) {
					add_parameter_from_stream(cur_key, value);
					value.str(std::string());
			}
			cur_key = argv[i];
		} else {
			if (!value.str().empty()) {
				value << " ";
			}
			value << argv[i];
		}
		i++;
	}
	add_parameter_from_stream(cur_key, value);
}

inline void Parser::add_parameter_from_stream(const char *key, const std::stringstream& value)
{
	parameters[std::string(key + 1)] = value.str();
}

std::ostream& operator<<(std::ostream& os, const Parser& parser)
{
	for (const auto& [key, value] : parser.parameters) {
		os << key << " : " << value << "\n";
	}
	return os;
}
