#include <iostream>
#include <argparse/argparse.hpp>
#include <PgmFileManager.hpp>
#include <GameOfLife.hpp>

void setup_parser(argparse::ArgumentParser& program)
{
	program.add_argument("-i")
		.help("initialize a grid")
		.default_value(false)
		.implicit_value(true);

	program.add_argument("-r")
		.help("run a simulation")
		.default_value(false)
		.implicit_value(true);

	program.add_argument("-k")
		.scan<'u', unsigned long>()
		//.default_value(std::string{"15"})
		.help("grid size");

	program.add_argument("-e")
		.scan<'u', unsigned char>()
		//.default_value(std::string{"1"})
		.help("evolution type (0 = ordered, 1 = static)");

	program.add_argument("-f")
		.default_value(std::string{"grid.pgm"})
		.help("input file name");

	program.add_argument("-n")
		.scan<'u', unsigned int>()
		//.default_value(15)
		.help("simulation steps");

	program.add_argument("-s")
		.scan<'u', unsigned int>()
		//.default_value(1)
		.help("snapshotting period");
}

int main(int argc, char **argv)
{
	int ret = EXIT_SUCCESS;
	argparse::ArgumentParser program{"game_of_life"};
	setup_parser(program);

	try {
		program.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << program;
		std::exit(1);
	}

	if (program["-i"] == true && program["-r"] == false) {
		auto filename = program.get<std::string>("-f");
		auto size = program.get<unsigned long>("-k");
		PgmFileManager pgm_manager{filename, size};
		pgm_manager.write();
	} else if (program["-r"] == true && program["-i"] == false) {
		auto filename = program.get<std::string>("-f");
		auto steps = program.get<unsigned int>("-n");
		auto evolution_strategy = program.get<unsigned char>("-e");
		auto snapshotting_period = program.get<unsigned int>("-s");
		PgmFileManager pgm_manager{filename};
		GameOfLife game{evolution_strategy, steps, snapshotting_period, pgm_manager};
		game.evolve();
	} else {
		std::cerr << "ERROR: Exactly one between -i and -r must be specified. Quitting." << std::endl;
		ret = EXIT_FAILURE;
	}

	return ret;
}
