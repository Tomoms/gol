#ifndef CELL_H
#define CELL_H

#include <vector>
#include <ostream>

#define CELL_ALIVE_VALUE	0
#define CELL_DEAD_VALUE		255

class Cell
{
private:
	unsigned long x_, y_;
	std::vector<std::reference_wrapper<Cell>> neighbors_;
	bool alive_;

public:
	Cell(unsigned long x, unsigned long y, bool alive);
	void add_to_neighbors(std::initializer_list<std::reference_wrapper<Cell>> cells);
	std::vector<std::reference_wrapper<Cell>>& get_neighbors();
	unsigned long get_x() const;
	unsigned long get_y() const;
	bool is_alive() const;
	bool is_alive_after_evolution() const;
	void set_alive(bool alive);
	friend std::ostream& operator<<(std::ostream& os, Cell const & cell);
};

#endif
