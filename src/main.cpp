#include <climits>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <argparse/argparse.hpp>
#include <boost/mpi.hpp>
#include <PgmUtils.hpp>
#include <GameOfLife.hpp>
#include <mpi.h>

#ifdef DEBUG
#define ALL_RANKS_PRINT(x) \
	std::cout << "Rank " << world.rank() << ": " << x << std::endl
#else
#define ALL_RANKS_PRINT(x) do {} while (0)
#endif

#ifdef DEBUG
#define ONE_RANK_PRINTS(r, x) \
	if (world.rank() == r) std::cout << "Rank " << r << ": " << x << std::endl
#else
#define ONE_RANK_PRINTS(r, x) do {} while (0)
#endif

namespace mpi = boost::mpi;

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

std::pair<ulong, ulong> compute_rank_chunk_bounds(ulong grid_size, mpi::communicator world)
{
	const auto cells = grid_size * grid_size;
	auto rank_elements = cells / world.size();
	auto regular_elements = rank_elements;
	const auto leftovers = cells % world.size();
	bool correction_factor = ulong(world.rank()) < leftovers;
	rank_elements += correction_factor;
	auto rank_offset = 0;
	for (auto i = 0; i < world.rank(); i++) {
		rank_offset += regular_elements + (ulong(i) < leftovers);
	}
	return { rank_elements, rank_offset };
}

int main(int argc, char **argv)
{
	mpi::environment env(argc, argv);
	mpi::communicator world;
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

	const auto filename = program.get<std::string>("-f");

	if (program["-i"] == true && program["-r"] == false) {
		ALL_RANKS_PRINT("enters image generation branch");
		const auto grid_size = program.get<unsigned long>("-k");
		auto [rank_elements, rank_offset] = compute_rank_chunk_bounds(grid_size, world);
		ALL_RANKS_PRINT("works on " << rank_elements << " elements");

		PGM_HOLDER rank_random_chunk = PgmUtils::generate_random_chunk(rank_elements);

		if (!world.rank()) {
			const SIZE_HOLDER dimensions{grid_size, grid_size};
			PgmUtils::write_header(filename, dimensions);
		}

		world.barrier();
		auto file_size = std::filesystem::file_size(filename);
		world.barrier();
		rank_offset += file_size;
		std::streampos rank_offset_streampos = static_cast<std::streampos>(rank_offset);
		ALL_RANKS_PRINT("offset " << rank_offset_streampos);
		PgmUtils::write_chunk_to_file(filename, rank_random_chunk, rank_offset_streampos, static_cast<MPI_Comm>(world));
	} else if (program["-i"] == false && program["-r"] == true) {
		ALL_RANKS_PRINT("enters running branch\nFirst step: read the image");
		ulong grid_size;
		uint header_length;

		// Rank 0 reads the header
		if (!world.rank()) {
			std::ifstream infile(filename.c_str());
			std::string line;
			std::getline(infile, line);
			std::istringstream iss(line);
			header_length = uint(line.size()) + 1; // account for new line
			std::string magic;
			iss >> magic >> grid_size;
			ALL_RANKS_PRINT(magic);
		}
		broadcast(world, grid_size, 0);
		broadcast(world, header_length, 0);

		ALL_RANKS_PRINT("detected size of " << grid_size);
		ALL_RANKS_PRINT("header length: " << header_length);

		auto [rank_elements, rank_offset] = compute_rank_chunk_bounds(grid_size, world);
		rank_offset += header_length;
		std::streampos rank_offset_streampos = static_cast<std::streampos>(rank_offset);
		ALL_RANKS_PRINT("offset " << rank_offset_streampos);
		PGM_HOLDER rank_chunk = PgmUtils::read_chunk_from_file(filename, rank_elements, rank_offset_streampos, static_cast<MPI_Comm>(world));
	} else {
		ONE_RANK_PRINTS(0, "invalid arguments, quitting.");
		ret = EXIT_FAILURE;
	}

	return ret;
}
