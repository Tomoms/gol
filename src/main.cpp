#include <iostream>
#include <argparse/argparse.hpp>
#include <boost/mpi.hpp>
#include <PgmUtils.hpp>
#include <GameOfLife.hpp>

namespace mpi = boost::mpi;

const PGM_HOLDER example_image{
	{
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
	}
};

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

std::string compute_checkpoint_filename(unsigned long step)
{
	std::string suffix{};
	if (step < 10) {
		suffix = "0000" + std::to_string(step);
	} else if (step < 100) {
		suffix = "000" + std::to_string(step);
	} else if (step < 1000) {
		suffix = "00" + std::to_string(step);
	} else if (step < 10000) {
		suffix = "0" + std::to_string(step);
	} else {
		suffix = std::to_string(step);
	}
	return "snapshot_" + suffix;
}

int main(int argc, char **argv)
{
	mpi::environment env(argc, argv);
	mpi::communicator world;
	int ret = EXIT_SUCCESS;
	bool run_to_generate_image = 0;
	if (world.rank() == 0) {
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
			run_to_generate_image = 1;
			broadcast(world, run_to_generate_image, 0);
			auto filename = program.get<std::string>("-f");
			auto size = program.get<unsigned long>("-k");
			// TODO implement generating a new image
		} else if (program["-r"] == true && program["-i"] == false) {
			broadcast(world, run_to_generate_image, 0);
			auto filename = program.get<std::string>("-f");
			auto steps = program.get<unsigned int>("-n");
			auto evolution_strategy = program.get<unsigned char>("-e");
			auto snapshotting_period = program.get<unsigned int>("-s");
			if (!snapshotting_period) {
				snapshotting_period = steps;
			}
			const SIZE_HOLDER grid_size = PgmUtils::read_size(filename);
			unsigned long elements = grid_size.first * grid_size.second;
			broadcast(world, elements, 0);
			const PGM_HOLDER image_data = PgmUtils::read_image_data(filename, grid_size);
			GameOfLife game{evolution_strategy, steps, snapshotting_period, grid_size, image_data};
			for (auto i = 1UL; i <= steps; i++) {
				game.evolve();
				if (i % snapshotting_period == 0) {
					std::string checkpoint_filename{compute_checkpoint_filename(i)};
					PgmUtils::write_image_data(checkpoint_filename, grid_size, PgmUtils::grid_to_image_data(game.get_grid()));
				}
			}
		} else {
			std::cerr << "ERROR: Exactly one between -i and -r must be specified. Quitting." << std::endl;
			ret = EXIT_FAILURE;
		}
	} else {
		broadcast(world, run_to_generate_image, 0);
		if (run_to_generate_image) {
			return ret;
		}
		unsigned long elements;
		broadcast(world, elements, 0);
	}

	return ret;
}
