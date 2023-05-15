#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <PgmFileManager.hpp>

class GameOfLife
{
private:
	unsigned char evolution_strategy_;
	unsigned int steps_;
	unsigned int snapshotting_period_;
	std::vector<bool> grid_;

public:
	GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, const PGM_HOLDER& image_data);
};

#endif
