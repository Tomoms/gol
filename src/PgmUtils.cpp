#include <limits>
#include <PgmUtils.hpp>

void PgmUtils::write_header(const std::string& filename, const SIZE_HOLDER& dimensions)
{
	std::ofstream outstream{filename.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << dimensions.first << " " << dimensions.second << " " << PGM_MAX_VALUE << std::endl;
}

void PgmUtils::write_chunk_to_file(const std::string& filename, const std::vector<char>& chunk,
								   std::streampos start_offset) {
    std::ofstream file(filename, std::ios::binary | std::ios::in | std::ios::out);
    file.seekp(start_offset);
    for (const char byte : chunk) {
        file.put(byte);
    }
    file.close();
}
