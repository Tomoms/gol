#ifndef PGMUTILS_H
#define PGMUTILS_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <Cell.hpp>

#define PGM_MAX_VALUE	255
#define PGM_HOLDER		std::vector<unsigned char>
#define SIZE_HOLDER		std::pair<unsigned long, unsigned long>

namespace PgmUtils {

SIZE_HOLDER read_size(std::string& filename);
PGM_HOLDER read_image_data(std::string& filename);
void write_image_data(std::string& filename, const SIZE_HOLDER& size, const PGM_HOLDER& image_data);
PGM_HOLDER grid_to_image_data(std::vector<Cell>& grid);

}

#endif
