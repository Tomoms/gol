#include <limits>
#include <PgmUtils.hpp>

SIZE_HOLDER PgmUtils::read_size(std::string& filename)
{
	unsigned long rows, cols;
	std::ifstream instream(filename.c_str(), std::ios_base::binary);
	std::string magic;
	instream >> magic >> cols >> rows;
	return { rows, cols };
}

PGM_HOLDER PgmUtils::read_image_data(std::string& filename, const SIZE_HOLDER& size)
{
	std::ifstream instream(filename.c_str(), std::ios_base::binary);
	instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	PGM_HOLDER image_data(size.first * size.second, 0);
	instream.read(reinterpret_cast<char*>(&image_data[0]), size.first * size.second);
	return image_data;
}

void PgmUtils::write_image_data(std::string& filename, const SIZE_HOLDER& size, const PGM_HOLDER& image_data)
{
	std::ofstream outstream{filename.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << size.second << " " << size.first << " " << PGM_MAX_VALUE << std::endl;
	outstream.write(reinterpret_cast<const char*>(&image_data[0]), size.first * size.second);
}

PGM_HOLDER PgmUtils::grid_to_image_data(std::vector<Cell>& grid)
{
	PGM_HOLDER image_data(grid.size(), 0);
	for (auto i = 0UL; i < grid.size(); i++) {
		if (grid[i].is_alive()) {
			image_data[i] = CELL_ALIVE_VALUE;
		} else {
			image_data[i] = CELL_DEAD_VALUE;
		}
	}
	return image_data;
}
