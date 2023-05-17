#include <iostream>
#include <GameOfLife.hpp>

#define REFERENCE_TO_NEIGHBOR(x, y) \
		std::ref(grid[coords_to_index(x, y)])

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
	populate_neighbors(grid_);
}

void GameOfLife::populate_neighbors(std::vector<Cell>& grid)
{
	for (auto& cell : grid) {
		unsigned long x = cell.get_x();
		unsigned long y = cell.get_y();
		unsigned long prev_row = x != 0 ? x - 1 : rows_ - 1;
		unsigned long prev_col = y != 0 ? y - 1 : cols_ - 1;
		unsigned long next_row = x != rows_ - 1 ? x + 1 : 0;
		unsigned long next_col = y != cols_ - 1 ? y + 1 : 0;
		cell.get_neighbors().clear();
		cell.add_to_neighbors(
		{
			REFERENCE_TO_NEIGHBOR(prev_row, prev_col),
			REFERENCE_TO_NEIGHBOR(prev_row, y),
			REFERENCE_TO_NEIGHBOR(prev_row, next_col),
			REFERENCE_TO_NEIGHBOR(x, prev_col),
			REFERENCE_TO_NEIGHBOR(x, next_col),
			REFERENCE_TO_NEIGHBOR(next_row, prev_col),
			REFERENCE_TO_NEIGHBOR(next_row, y),
			REFERENCE_TO_NEIGHBOR(next_row, next_col),
		}
		);
	}
}

void GameOfLife::evolve()
{
	if (evolution_strategy_) { // static
		std::vector<Cell> new_grid = grid_;
		populate_neighbors(new_grid);
		for (auto i = 0UL; i < grid_.size(); i++) {
			bool new_status = grid_[i].becomes_alive();
			new_grid[i].set_alive(new_status);
		}
	} else { // ordered
	}
}

inline unsigned long GameOfLife::coords_to_index(unsigned long x, unsigned long y)
{
	return x * cols_ + y;
}
