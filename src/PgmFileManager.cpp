#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <PgmFileManager.hpp>

PgmFileManager::PgmFileManager(std::string& filename, unsigned long size):
filename_{filename}, rows_{size}, cols_{size}
{}

PgmFileManager::PgmFileManager(std::string& filename):
filename_{filename}
{
	{
		std::ifstream instream(filename.c_str(), std::ios_base::binary);
		std::string magic;
		instream >> magic >> rows_ >> cols_;
	}

	image_data_.resize(rows_ * cols_, 0);

	{
		std::ifstream instream(filename.c_str(), std::ios_base::binary);
		instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		instream.read(reinterpret_cast<char*>(&image_data_[0]), rows_ * cols_);
	}
}

void PgmFileManager::write(const PGM_HOLDER& image_data)
{
	std::ofstream outstream{filename_.c_str(), std::ios_base::binary | std::ios_base::trunc};
	outstream << "P5 " << cols_ << " " << rows_ << " " << PGM_MAX_VALUE << std::endl;
	outstream.write(reinterpret_cast<const char*>(&image_data[0]), rows_ * cols_);
}

unsigned long PgmFileManager::get_rows()
{
	return rows_;
}

unsigned long PgmFileManager::get_cols()
{
	return cols_;
}

PGM_HOLDER& PgmFileManager::get_image_data()
{
	return image_data_;
}

PGM_HOLDER PgmFileManager::grid_to_image_data(std::vector<Cell>& grid)
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
