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
	const auto base_rows_per_rank = grid_size / world.size();
	auto rank_rows = base_rows_per_rank;
	const auto leftovers = grid_size % world.size();
	if (ulong(world.rank()) < leftovers)
		rank_rows++;
	auto rank_offset = 0UL;
	for (auto i = 0; i < world.rank(); i++) {
		bool correction_factor = ulong(i) < leftovers;
		rank_offset += (base_rows_per_rank + correction_factor) * grid_size;
	}
	return { rank_rows, rank_offset };
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
		const auto grid_size = program.get<unsigned long>("-k");
		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(grid_size, world);
		ALL_RANKS_PRINT("works on " << rank_rows * grid_size << " elements");

		PGM_HOLDER rank_random_chunk = PgmUtils::generate_random_chunk(rank_rows * grid_size);

		if (!world.rank()) {
			const SIZE_HOLDER dimensions{grid_size, grid_size};
			PgmUtils::write_header(filename, dimensions);
		}

		world.barrier();
		auto file_size = std::filesystem::file_size(filename);
		world.barrier();
		auto rank_file_offset = rank_offset + file_size;
		std::streampos rank_file_offset_streampos = static_cast<std::streampos>(rank_file_offset);
		ALL_RANKS_PRINT("offset " << rank_file_offset_streampos);
		PgmUtils::write_chunk_to_file(filename, rank_random_chunk, rank_file_offset_streampos, 0, static_cast<MPI_Comm>(world));
	} else if (program["-i"] == false && program["-r"] == true) {
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
		}
		broadcast(world, grid_size, 0);
		broadcast(world, header_length, 0);

		ALL_RANKS_PRINT("detected size of " << grid_size);
		ALL_RANKS_PRINT("header length: " << header_length);

		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(grid_size, world);
		ALL_RANKS_PRINT("will work on " << rank_rows << " rows, i.e. " << rank_rows * grid_size << " cells");
		auto rank_file_offset = rank_offset + header_length;
		std::streampos rank_file_offset_streampos = static_cast<std::streampos>(rank_file_offset);
		ALL_RANKS_PRINT("offset " << rank_file_offset_streampos);
		PGM_HOLDER rank_chunk = PgmUtils::read_chunk_from_file(filename, rank_rows * grid_size, rank_file_offset_streampos, grid_size, static_cast<MPI_Comm>(world));

		int prev_rank = world.rank() - 1 >= 0 ? world.rank() - 1 : world.size() - 1;
		int next_rank = world.rank() + 1 >= world.size() ? 0 : world.rank() + 1;
		if (world.rank()) {
			world.send(prev_rank, 1, rank_chunk.data() + grid_size, grid_size);
			world.send(next_rank, 2, rank_chunk.data() + rank_rows * grid_size, grid_size);
			world.recv(prev_rank, 2, rank_chunk.data(), grid_size);
			world.recv(next_rank, 1, rank_chunk.data() + (rank_rows + 1) * grid_size, grid_size);
		} else {
			world.recv(prev_rank, 2, rank_chunk.data(), grid_size);
			world.recv(next_rank, 1, rank_chunk.data() + (rank_rows + 1) * grid_size, grid_size);
			world.send(prev_rank, 1, rank_chunk.data() + grid_size, grid_size);
			world.send(next_rank, 2, rank_chunk.data() + rank_rows * grid_size, grid_size);
		}
	} else {
		ONE_RANK_PRINTS(0, "invalid arguments, quitting.");
		ret = EXIT_FAILURE;
	}

	return ret;
}
