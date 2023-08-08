#ifndef PGMUTILS_H
#define PGMUTILS_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <Cell.hpp>
#include <mpi.h>

#define PGM_MAX_VALUE	255
#define PGM_HOLDER		std::vector<unsigned char>
#define SIZE_HOLDER		std::pair<unsigned long, unsigned long> // width (number of columns) and height (number of rows)

namespace PgmUtils {

	void write_header(const std::string& filename, const SIZE_HOLDER& dimensions);
	void write_chunk_to_file(const std::string& filename, const PGM_HOLDER& chunk,
							 const std::streampos start_offset, MPI_Comm comm);
}

#endif
