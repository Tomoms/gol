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
#include <omp.h>

#ifdef DEBUG
#define ALL_RANKS_PRINT(x) \
	std::cout << "Rank " << world.rank() << ": " << x << std::endl
#define ONE_RANK_PRINTS(r, x) \
	if (world.rank() == r) std::cout << "Rank " << r << ": " << x << std::endl
#else
#define ALL_RANKS_PRINT(x) do {} while (0)
#define ONE_RANK_PRINTS(r, x) do {} while (0)
#endif

#define FIRST_ROW_OF_SENDING_RANK	1
#define LAST_ROW_OF_SENDING_RANK	2

#define CELL_ALIVE	255
#define CELL_DEAD	0
#define IS_CELL_ALIVE(index) rank_chunk[index] == CELL_ALIVE
#define IS_TOP_NEIGHBOR_ALIVE(index) IS_CELL_ALIVE(index - grid_size)
#define IS_BOTTOM_NEIGHBOR_ALIVE(index) IS_CELL_ALIVE(index + grid_size)

#define SEND_LAST_ROW \
	world.send(next_rank, LAST_ROW_OF_SENDING_RANK, rank_chunk.data() + rank_rows * grid_size, grid_size);
#define SEND_FIRST_ROW \
	world.send(prev_rank, FIRST_ROW_OF_SENDING_RANK, rank_chunk.data() + grid_size, grid_size);
#define RECEIVE_LAST_ROW \
	world.recv(prev_rank, LAST_ROW_OF_SENDING_RANK, rank_chunk.data(), grid_size);
#define RECEIVE_FIRST_ROW \
	world.recv(next_rank, FIRST_ROW_OF_SENDING_RANK, rank_chunk.data() + (rank_rows + 1) * grid_size, grid_size);

namespace mpi = boost::mpi;
namespace mt  = mpi::threading;

int prev_rank, next_rank;
ulong grid_size;

inline bool is_top_left_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if (index % grid_size == 0) { // left edge of the matrix
		alive = IS_CELL_ALIVE(index - 1);
	} else {
		alive = IS_CELL_ALIVE(index - grid_size - 1);
	}
	return alive;
}

inline bool is_top_right_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if ((index + 1) % grid_size == 0) { // right edge of the matrix
		alive = IS_CELL_ALIVE(index + 1 - 2 * grid_size);
	} else {
		alive = IS_CELL_ALIVE(index - grid_size + 1);
	}
	return alive;
}

inline bool is_left_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if (index % grid_size == 0) { // left edge of the matrix
		alive = IS_CELL_ALIVE(index + grid_size - 1);
	} else {
		alive = IS_CELL_ALIVE(index - 1);
	}
	return alive;
}

inline bool is_right_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if ((index + 1) % grid_size == 0) { // right edge of the matrix
		alive = IS_CELL_ALIVE(index - grid_size + 1);
	} else {
		alive = IS_CELL_ALIVE(index + 1);
	}
	return alive;
}

inline bool is_bottom_left_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if (index % grid_size == 0) { // left edge of the matrix
		alive = IS_CELL_ALIVE(index - 1);
	} else {
		alive = IS_CELL_ALIVE(index + 2 * grid_size - 1);
	}
	return alive;
}

inline bool is_bottom_right_neighbor_alive(PGM_HOLDER& rank_chunk, ulong index)
{
	bool alive = 0;
	if ((index + 1) % grid_size == 0) { // right edge of the matrix
		alive = IS_CELL_ALIVE(index + 1);
	} else {
		alive = IS_CELL_ALIVE(index + grid_size + 1);
	}
	return alive;
}

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

std::pair<ulong, ulong> compute_rank_chunk_bounds(mpi::communicator world)
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

inline char count_alive_neighbors(PGM_HOLDER& rank_chunk, ulong j)
{
	return is_bottom_left_neighbor_alive(rank_chunk, j) + IS_BOTTOM_NEIGHBOR_ALIVE(j) +
		is_bottom_right_neighbor_alive(rank_chunk, j) + is_left_neighbor_alive(rank_chunk, j) +
		is_right_neighbor_alive(rank_chunk, j) + is_top_left_neighbor_alive(rank_chunk, j) +
		IS_TOP_NEIGHBOR_ALIVE(j) + is_top_right_neighbor_alive(rank_chunk, j);
}

PGM_HOLDER evolve_static(PGM_HOLDER& rank_chunk, mpi::communicator world, int i)
{
	const ulong rank_rows = (rank_chunk.size() / grid_size) - 2; // minus 2 halo rows
	PGM_HOLDER next_step_chunk((rank_rows + 2) * grid_size);
	if (world.rank()) {
		SEND_FIRST_ROW;
		SEND_LAST_ROW;
		RECEIVE_LAST_ROW;
		RECEIVE_FIRST_ROW;
	} else {
		RECEIVE_LAST_ROW;
		RECEIVE_FIRST_ROW;
		SEND_FIRST_ROW;
		SEND_LAST_ROW;
	}

	int rank = world.rank();
	for (auto startpos = grid_size; startpos < rank_chunk.size() - grid_size; startpos += grid_size) {
		#pragma omp task shared(rank_chunk, next_step_chunk, grid_size)
		{
			if (i==14)
				std::cout << "rank " << rank << " starting at " << startpos << ", executed by " << omp_get_thread_num() << std::endl;
			for (auto j = startpos; j < startpos + grid_size; j++) {
				char alive_neighbors = count_alive_neighbors(rank_chunk, j);
				if (alive_neighbors == 3) {
					next_step_chunk[j] = CELL_ALIVE;
				} else if (alive_neighbors == 2) {
					next_step_chunk[j] = rank_chunk[j];
				} else {
					next_step_chunk[j] = CELL_DEAD;
				}
			}
		}

	}
	return next_step_chunk;
}

inline void update_cell_ordered(PGM_HOLDER& rank_chunk, ulong j)
{
	char alive_neighbors = count_alive_neighbors(rank_chunk, j);
	if (alive_neighbors == 3) {
		rank_chunk[j] = CELL_ALIVE;
	} else if (alive_neighbors != 2) {
		rank_chunk[j] = CELL_DEAD;
	}
}

PGM_HOLDER evolve_ordered(PGM_HOLDER& rank_chunk, mpi::communicator world, int i)
{
	const ulong rank_rows = (rank_chunk.size() / grid_size) - 2;
	if (world.rank() == 0) {
		RECEIVE_LAST_ROW;
		RECEIVE_FIRST_ROW;
		for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
			update_cell_ordered(rank_chunk, j);
		}
		SEND_FIRST_ROW;
		SEND_LAST_ROW;
	} else if (world.rank() == world.size() - 1) {
		SEND_LAST_ROW;
		SEND_FIRST_ROW;
		RECEIVE_FIRST_ROW;
		RECEIVE_LAST_ROW;
		for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
			update_cell_ordered(rank_chunk, j);
		}
	} else {
		SEND_FIRST_ROW;
		RECEIVE_FIRST_ROW;
		RECEIVE_LAST_ROW;
		for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
			update_cell_ordered(rank_chunk, j);
		}
		SEND_LAST_ROW;
	}
	return rank_chunk;
}

void save_snapshot(PGM_HOLDER& rank_chunk, int i, std::streampos rank_file_offset_streampos, mpi::communicator world)
{
	const auto checkpoint_filename = compute_checkpoint_filename(i);
	if (!world.rank()) {
		const SIZE_HOLDER dimensions{grid_size, grid_size};
		PgmUtils::write_header(checkpoint_filename, dimensions);
	}
	PgmUtils::write_chunk_to_file(checkpoint_filename, rank_chunk, rank_file_offset_streampos, grid_size, static_cast<MPI_Comm>(world));
}

int main(int argc, char **argv)
{
	mpi::environment env(argc, argv, mt::funneled);
	if (env.thread_level() < mt::funneled) {
		env.abort(-1);
	}
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
		grid_size = program.get<unsigned long>("-k");
		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(world);
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

		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(world);
		ALL_RANKS_PRINT("will work on " << rank_rows << " rows, i.e. " << rank_rows * grid_size << " cells");
		auto rank_file_offset = rank_offset + header_length;
		std::streampos rank_file_offset_streampos = static_cast<std::streampos>(rank_file_offset);
		ALL_RANKS_PRINT("offset " << rank_file_offset_streampos);
		PGM_HOLDER rank_chunk = PgmUtils::read_chunk_from_file(filename, rank_rows * grid_size, rank_file_offset_streampos, grid_size, static_cast<MPI_Comm>(world));

		prev_rank = world.rank() - 1 >= 0 ? world.rank() - 1 : world.size() - 1;
		next_rank = world.rank() + 1 >= world.size() ? 0 : world.rank() + 1;

		const auto simulation_steps = program.get<unsigned int>("-n");
		const auto snapshotting_period = program.get<unsigned int>("-s");
		const auto evolution_type = program.get<unsigned char>("-e");

		PGM_HOLDER (*evolver)(PGM_HOLDER&, mpi::communicator, int);
		if (evolution_type == 1) {
			evolver = evolve_static;
		} else if (evolution_type == 0) {
			evolver = evolve_ordered;
		} else {
			ONE_RANK_PRINTS(0, "Unknown evolution type. Quitting.");
			ret = EXIT_FAILURE;
			return ret;
		}

		#pragma omp parallel
		{
		#pragma omp master
			{
				for (uint i = 1; i <= simulation_steps; i++) {
					rank_chunk = evolver(rank_chunk, world, i);
					if (snapshotting_period) {
						if (i % snapshotting_period == 0) {
							save_snapshot(rank_chunk, i, rank_file_offset_streampos, world);
						}
					} else {
						if (i == simulation_steps) {
							save_snapshot(rank_chunk, i, rank_file_offset_streampos, world);
						}
					}
				}
			}
		}
	} else {
		ONE_RANK_PRINTS(0, "invalid arguments, quitting.");
		ret = EXIT_FAILURE;
	}

	return ret;
}
