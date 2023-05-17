#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <PgmFileManager.hpp>
#include <Cell.hpp>

#define CELL_ALIVE_VALUE	0
#define CELL_DEAD_VALUE		255

class GameOfLife
{
private:
	unsigned char evolution_strategy_;
	unsigned int steps_;
	unsigned int snapshotting_period_;
	std::vector<Cell> grid_;
	unsigned long rows_;
	unsigned long cols_;

	unsigned long coords_to_index(unsigned long x, unsigned long y);
	void populate_neighbors(std::vector<Cell>& grid);

public:
	GameOfLife(unsigned char evolution_strategy, unsigned int steps, unsigned int snapshotting_period, PgmFileManager& pgm_manager);
	void evolve();
	void print_grid(std::vector<Cell>& grid) const;
	std::vector<Cell>& get_grid(void);
};

#endif
