#include <limits>
#include <PgmUtils.hpp>

void PgmUtils::write_header(const std::string& filename, const SIZE_HOLDER& dimensions)
{
	std::ofstream outstream{filename.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << dimensions.first << " " << dimensions.second << " " << PGM_MAX_VALUE << std::endl;
}
