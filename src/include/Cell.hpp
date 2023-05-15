#ifndef CELL_H
#define CELL_H

#include <vector>

class Cell
{
private:
	unsigned long x_, y_;
	std::vector<std::reference_wrapper<Cell>> neighbors_;
	bool alive_;

public:
	Cell(unsigned long x, unsigned long y, bool alive);
	void add_to_neighbors(std::initializer_list<std::reference_wrapper<Cell>> cells);
	unsigned long get_x();
	unsigned long get_y();
};

#endif
