#include <limits>
#include <PgmUtils.hpp>
#include <iostream>

void PgmUtils::write_header(const std::string& filename, const SIZE_HOLDER& dimensions)
{
	std::ofstream outstream{filename.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << dimensions.first << " " << dimensions.second << " " << PGM_MAX_VALUE << std::endl;
}

void PgmUtils::write_chunk_to_file(const std::string& filename, const PGM_HOLDER& chunk,
									const std::streampos start_offset, const ulong leading_halo_length,
									MPI_Comm comm)
{
	MPI_File file;
	MPI_File_open(comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	std::size_t size = chunk.size() - 2 * leading_halo_length;
	MPI_Offset offset = static_cast<MPI_Offset>(start_offset);
	MPI_File_write_at_all(file, offset, chunk.data() + leading_halo_length, size, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&file);
}

PGM_HOLDER PgmUtils::read_chunk_from_file(const std::string& filename, const ulong chunk_length,
									const std::streampos start_offset, const ulong leading_halo_length,
									MPI_Comm comm)
{
	MPI_File file;
	MPI_File_open(comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	PGM_HOLDER chunk(chunk_length + 2 * leading_halo_length);
	MPI_Offset offset = static_cast<MPI_Offset>(start_offset);
	MPI_File_read_at_all(file, offset, chunk.data() + leading_halo_length, chunk_length, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&file);
	return chunk;
}

PGM_HOLDER PgmUtils::generate_random_chunk(unsigned long size)
{
	PGM_HOLDER chunk(size);
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_int_distribution<int> distribution(0, 1);
	for (auto i = 0UL; i < size; i++) {
		int value = distribution(generator);
		chunk[i] = (value == 0) ? 0x00 : 0xFF;
	}
	return chunk;
}
