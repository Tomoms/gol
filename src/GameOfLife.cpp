#include <iostream>
#include <GameOfLife.hpp>

GameOfLife::GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, PGM_HOLDER& image_data):
evolution_strategy{evolution_strategy},
steps{steps},
snapshotting_period{snapshotting_period},
grid(image_data.size() * image_data[0].size(), 0)
{
	for (auto i = 0UL; i < grid.size(); i++) {
		auto row = i / image_data[0].size();
		auto col = i % image_data[0].size();
		grid[i] = image_data[row][col];
	}
}
