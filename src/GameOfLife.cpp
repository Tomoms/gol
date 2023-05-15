#include <iostream>
#include <GameOfLife.hpp>

GameOfLife::GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, PgmFileManager& pgm_manager):
evolution_strategy_{evolution_strategy},
steps_{steps},
snapshotting_period_{snapshotting_period},
rows_{pgm_manager.get_rows()},
cols_{pgm_manager.get_cols()}
{
	for (auto i = 0UL; i < rows_; i++) {
		for (auto j = 0UL; j < cols_; j++) {
			grid_.emplace_back(Cell{i, j, pgm_manager.get_image_data()[i][j] == CELL_ALIVE_VALUE});
		}
	}
	populate_neighbors();
}

void GameOfLife::populate_neighbors()
{
	for (auto& cell : grid_) {
		unsigned long x = cell.get_x();
		unsigned long y = cell.get_y();
		unsigned long prev_row = x != 0 ? x - 1 : rows_ - 1;
		unsigned long prev_col = y != 0 ? y - 1 : cols_ - 1;
		unsigned long next_row = x != rows_ - 1 ? x + 1 : 0;
		unsigned long next_col = y != cols_ - 1 ? y + 1 : 0;
		cell.add_to_neighbors(
		{
			std::ref(grid_[coords_to_index(prev_row, prev_col)]),
			std::ref(grid_[coords_to_index(prev_row, y)]),
			std::ref(grid_[coords_to_index(prev_row, next_col)]),
			std::ref(grid_[coords_to_index(x, prev_col)]),
			std::ref(grid_[coords_to_index(x, next_col)]),
			std::ref(grid_[coords_to_index(next_row, prev_col)]),
			std::ref(grid_[coords_to_index(next_row, y)]),
			std::ref(grid_[coords_to_index(next_row, next_col)]),
		}
		);
	}
}

inline unsigned long GameOfLife::coords_to_index(unsigned long x, unsigned long y)
{
	return x * cols_ + y;
}
