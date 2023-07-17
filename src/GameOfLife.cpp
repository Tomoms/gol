#include <iostream>
#include <GameOfLife.hpp>

#define REFERENCE_TO_CELL(x, y) \
		std::ref(grid[coords_to_index(x, y)])

GameOfLife::GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, const SIZE_HOLDER& size, const PGM_HOLDER& image_data):
evolution_strategy_{evolution_strategy},
steps_{steps},
snapshotting_period_{snapshotting_period},
rows_{size.first},
cols_{size.second}
{
	for (auto i = 0UL; i < rows_; i++) {
		for (auto j = 0UL; j < cols_; j++) {
			grid_.emplace_back(Cell{i, j, image_data[i * cols_ + j] == CELL_ALIVE_VALUE});
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
			REFERENCE_TO_CELL(prev_row, prev_col),
			REFERENCE_TO_CELL(prev_row, y),
			REFERENCE_TO_CELL(prev_row, next_col),
			REFERENCE_TO_CELL(x, prev_col),
			REFERENCE_TO_CELL(x, next_col),
			REFERENCE_TO_CELL(next_row, prev_col),
			REFERENCE_TO_CELL(next_row, y),
			REFERENCE_TO_CELL(next_row, next_col),
		}
		);
	}
}

void GameOfLife::evolve()
{
	if (evolution_strategy_) { // static
		std::vector<Cell> new_grid = grid_;
		populate_neighbors(new_grid);
#pragma omp parallel for
		for (auto i = 0UL; i < grid_.size(); i++) {
			bool new_status = grid_[i].is_alive_after_evolution();
			new_grid[i].set_alive(new_status);
		}
		grid_ = std::move(new_grid);
	} else { // ordered
	}
}

void GameOfLife::print_grid() const
{
	for (auto i = 0UL; i < rows_; i++) {
		for (auto j = 0UL; j < cols_; j++) {
			std::cout << (grid_[i * cols_ + j].is_alive() ? "1" : "0") << " ";
		}
		std::cout << "\n";
	}
}

std::vector<Cell>& GameOfLife::get_grid()
{
	return grid_;
}

inline unsigned long GameOfLife::coords_to_index(unsigned long x, unsigned long y)
{
	return x * cols_ + y;
}
