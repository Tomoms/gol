#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <PgmFileManager.hpp>

class GameOfLife
{
private:
	unsigned char evolution_strategy;
	unsigned int steps;
	unsigned int snapshotting_period;
	std::vector<bool> grid;

public:
	GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, PGM_HOLDER& image_data);
};

#endif
