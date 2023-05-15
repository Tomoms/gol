#include <iostream>
#include <GameOfLife.hpp>

GameOfLife::GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, const PGM_HOLDER& image_data):
evolution_strategy_{evolution_strategy},
steps_{steps},
snapshotting_period_{snapshotting_period},
grid_(image_data.size() * image_data[0].size(), 0)
{
	for (auto i = 0UL; i < grid_.size(); i++) {
		auto row = i / image_data[0].size();
		auto col = i % image_data[0].size();
		grid_[i] = image_data[row][col];
	}
}
