#include <climits>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <argparse/argparse.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/timer.hpp>
#include <PgmUtils.hpp>
#include <mpi.h>
#include <omp.h>
#include <mimalloc.h>
#include <mimalloc-new-delete.h>

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
#define IS_CELL_ALIVE(index) (rank_chunk[index] == CELL_ALIVE)

#define SEND_LAST_ROW \
	world.isend(next_rank, LAST_ROW_OF_SENDING_RANK, rank_chunk.data() + rank_rows * grid_size, grid_size);
#define SEND_FIRST_ROW \
	world.isend(prev_rank, FIRST_ROW_OF_SENDING_RANK, rank_chunk.data() + grid_size, grid_size);
#define RECEIVE_TOP_HALO \
	world.recv(prev_rank, LAST_ROW_OF_SENDING_RANK, rank_chunk.data(), grid_size);
#define RECEIVE_BOTTOM_HALO \
	world.recv(next_rank, FIRST_ROW_OF_SENDING_RANK, rank_chunk.data() + (rank_rows + 1) * grid_size, grid_size);

namespace mpi = boost::mpi;
namespace mt  = mpi::threading;

int prev_rank, next_rank;
ulong grid_size;
uint nthreads;

inline __attribute__((always_inline)) unsigned char check_left_side(PGM_HOLDER& rank_chunk, ulong index)
{
	return (index % grid_size == 0) ? IS_CELL_ALIVE(index - 1) + IS_CELL_ALIVE(index + grid_size - 1) + IS_CELL_ALIVE(index - 1 + 2 * grid_size)
		: IS_CELL_ALIVE(index - grid_size - 1) + IS_CELL_ALIVE(index - 1) + IS_CELL_ALIVE(index + grid_size - 1);
}

inline __attribute__((always_inline)) unsigned char check_right_side(PGM_HOLDER& rank_chunk, ulong index)
{
	return ((index + 1) % grid_size == 0) ? IS_CELL_ALIVE(index + 1 - 2 * grid_size) + IS_CELL_ALIVE(index - grid_size + 1) + IS_CELL_ALIVE(index + 1)
		: IS_CELL_ALIVE(index - grid_size + 1) + IS_CELL_ALIVE(index + 1) + IS_CELL_ALIVE(index + grid_size + 1);
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
		.help("grid size");

	program.add_argument("-e")
		.scan<'u', unsigned char>()
		.help("evolution type (0 = ordered, 1 = static)");

	program.add_argument("-f")
		.default_value(std::string{"grid.pgm"})
		.help("input file name");

	program.add_argument("-n")
		.scan<'u', unsigned int>()
		.help("simulation steps");

	program.add_argument("-s")
		.scan<'u', unsigned int>()
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
	return "pics/snapshot_" + suffix;
}

std::pair<ulong, ulong> compute_rank_chunk_bounds(mpi::communicator world)
{
	if (world.size() == 1) {
		return { grid_size, 0 };
	}
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

inline __attribute__((always_inline)) char count_alive_neighbors(PGM_HOLDER& rank_chunk, ulong j)
{
	return check_left_side(rank_chunk, j) + check_right_side(rank_chunk, j) + IS_CELL_ALIVE(j + grid_size) + IS_CELL_ALIVE(j - grid_size);
}

void evolve_static(PGM_HOLDER& rank_chunk, PGM_HOLDER& next_step_chunk, mpi::communicator world)
{
	const ulong rank_rows = (rank_chunk.size() / grid_size) - 2; // minus 2 halo rows
	if (world.size() != 1) {
		if (world.rank()) {
			SEND_FIRST_ROW;
			SEND_LAST_ROW;
			RECEIVE_TOP_HALO;
			RECEIVE_BOTTOM_HALO;
		} else {
			RECEIVE_TOP_HALO;
			RECEIVE_BOTTOM_HALO;
			SEND_FIRST_ROW;
			SEND_LAST_ROW;
		}
	}
	const auto start = grid_size, end = (rank_rows + 1) * grid_size;
	const auto elements = end - start;
	const ulong elements_per_thread = elements / nthreads;
#pragma omp taskloop shared(rank_chunk, next_step_chunk) grainsize(1)
	for (auto j = start; j < nthreads * elements_per_thread; j += elements_per_thread) {
		auto i = j;
		for (; i < j + elements_per_thread; i += 4) {
			const auto next_1 = i+1, next_2 = i+2, next_3 = i+3;
			char alive_neighbors = count_alive_neighbors(rank_chunk, i);
			char alive_neighbors_1 = count_alive_neighbors(rank_chunk, next_1);
			char alive_neighbors_2 = count_alive_neighbors(rank_chunk, next_2);
			char alive_neighbors_3 = count_alive_neighbors(rank_chunk, next_3);
			next_step_chunk[i] = (alive_neighbors == 3 || alive_neighbors == 2) ? CELL_ALIVE : CELL_DEAD;
			next_step_chunk[next_1] = (alive_neighbors_1 == 3 || alive_neighbors_1 == 2) ? CELL_ALIVE : CELL_DEAD;
			next_step_chunk[next_2] = (alive_neighbors_2 == 3 || alive_neighbors_2 == 2) ? CELL_ALIVE : CELL_DEAD;
			next_step_chunk[next_3] = (alive_neighbors_3 == 3 || alive_neighbors_3 == 2) ? CELL_ALIVE : CELL_DEAD;
		}
		if (i != j + elements_per_thread) {
			i = i - 4 + 1;
			for (; i < j + elements_per_thread; i++) {
				char alive_neighbors = count_alive_neighbors(rank_chunk, i);
				next_step_chunk[i] = (alive_neighbors == 3 || alive_neighbors == 2) ? CELL_ALIVE : CELL_DEAD;
			}
		}
	}
	for (auto i = start + nthreads * elements_per_thread; i < end; i++) {
		char alive_neighbors = count_alive_neighbors(rank_chunk, i);
		next_step_chunk[i] = (alive_neighbors == 3 || alive_neighbors == 2) ? CELL_ALIVE : CELL_DEAD;
	}
}

inline __attribute__((always_inline)) void update_cell_ordered(PGM_HOLDER& rank_chunk, ulong j)
{
	char alive_neighbors = count_alive_neighbors(rank_chunk, j);
	rank_chunk[j] = (alive_neighbors == 3 || alive_neighbors == 2) ? CELL_ALIVE : CELL_DEAD;
}

void evolve_ordered(PGM_HOLDER& rank_chunk, PGM_HOLDER& unused, mpi::communicator world)
{
	const ulong rank_rows = (rank_chunk.size() / grid_size) - 2;
	if (world.size() != 1) {
		if (world.rank() == 0) {
			RECEIVE_TOP_HALO;
			RECEIVE_BOTTOM_HALO;
			for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
				update_cell_ordered(rank_chunk, j);
			}
			SEND_FIRST_ROW;
			SEND_LAST_ROW;
		} else if (world.rank() == world.size() - 1) {
			SEND_LAST_ROW;
			SEND_FIRST_ROW;
			RECEIVE_BOTTOM_HALO;
			RECEIVE_TOP_HALO;
			for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
				update_cell_ordered(rank_chunk, j);
			}
		} else {
			SEND_FIRST_ROW;
			RECEIVE_BOTTOM_HALO;
			RECEIVE_TOP_HALO;
			for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
				update_cell_ordered(rank_chunk, j);
			}
			SEND_LAST_ROW;
		}
	} else {
		for (auto j = grid_size; j < (rank_rows + 1) * grid_size ; j++) {
			update_cell_ordered(rank_chunk, j);
		}
	}
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
	const auto ranks = world.size();

	if (program["-i"] == true && program["-r"] == false) {
		grid_size = program.get<unsigned long>("-k");
		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(world);

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

		if (ranks != 1) {
			broadcast(world, grid_size, 0);
			broadcast(world, header_length, 0);
			prev_rank = world.rank() - 1 >= 0 ? world.rank() - 1 : world.size() - 1;
			next_rank = world.rank() + 1 >= world.size() ? 0 : world.rank() + 1;
		}

		const auto simulation_steps = program.get<unsigned int>("-n");
		const auto snapshotting_period = program.get<unsigned int>("-s");
		const auto evolution_type = program.get<unsigned char>("-e");

		void (*evolver)(PGM_HOLDER&, PGM_HOLDER&, mpi::communicator);
		if (evolution_type == 1) {
			evolver = evolve_static;
		} else if (evolution_type == 0) {
			evolver = evolve_ordered;
		} else {
			ONE_RANK_PRINTS(0, "Unknown evolution type. Quitting.");
			ret = EXIT_FAILURE;
			return ret;
		}

		auto [rank_rows, rank_offset] = compute_rank_chunk_bounds(world);
		auto rank_file_offset = rank_offset + header_length;
		std::streampos rank_file_offset_streampos = static_cast<std::streampos>(rank_file_offset);
		mpi::timer timer;
		PGM_HOLDER rank_chunk = PgmUtils::read_chunk_from_file(filename, rank_rows * grid_size, rank_file_offset_streampos, grid_size, static_cast<MPI_Comm>(world));
		PGM_HOLDER next_step_chunk(rank_chunk.size());

#pragma omp parallel
{
#pragma omp master
{
		nthreads = omp_get_num_threads();
		for (uint i = 1; i <= simulation_steps; i++) {
			evolver(rank_chunk, next_step_chunk, world);
			rank_chunk.swap(next_step_chunk);
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
		double elapsed = timer.elapsed();
		double avg = mpi::all_reduce(world, elapsed, std::plus<double>());
		avg = avg / world.size();
		if (!world.rank()){
			std::cout << grid_size << "," << world.size() << "," << nthreads << "," << avg << std::endl;
		}
	} else {
		ONE_RANK_PRINTS(0, "invalid arguments, quitting.");
		ret = EXIT_FAILURE;
	}

	return ret;
}
