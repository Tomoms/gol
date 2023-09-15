#ifndef PGMUTILS_H
#define PGMUTILS_H

#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <mpi.h>
#include <mimalloc.h>

#define PGM_MAX_VALUE	255
#define PGM_HOLDER		std::vector<unsigned char, mi_stl_allocator<unsigned char>>
#define SIZE_HOLDER		std::pair<unsigned long, unsigned long> // width (number of columns) and height (number of rows)

namespace PgmUtils {

	void write_header(const std::string& filename, const SIZE_HOLDER& dimensions);
	void write_chunk_to_file(const std::string& filename, const PGM_HOLDER& chunk,
							const std::streampos start_offset, const ulong leading_halo_length,
							MPI_Comm comm);
	PGM_HOLDER read_chunk_from_file(const std::string& filename, const ulong chunk_length,
									const std::streampos start_offset, const ulong leading_halo_length,
									MPI_Comm comm);
	PGM_HOLDER generate_random_chunk(unsigned long size);
}

#endif
