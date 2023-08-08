#include <limits>
#include <PgmUtils.hpp>

void PgmUtils::write_header(const std::string& filename, const SIZE_HOLDER& dimensions)
{
	std::ofstream outstream{filename.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << dimensions.first << " " << dimensions.second << " " << PGM_MAX_VALUE << std::endl;
}

void PgmUtils::write_chunk_to_file(const std::string& filename, const PGM_HOLDER& chunk,
									const std::streampos start_offset, MPI_Comm comm)
{
	MPI_File file;
	MPI_File_open(comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	std::size_t size = chunk.size();
	MPI_Offset offset = static_cast<MPI_Offset>(start_offset);
	MPI_File_write_at_all(file, offset, chunk.data(), size, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&file);
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
