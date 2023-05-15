#include <Cell.hpp>

Cell::Cell(unsigned long x, unsigned long y, bool alive):
x_{x},
y_{y},
neighbors_{},
alive_{alive}
{}

void Cell::add_to_neighbors(std::initializer_list<std::reference_wrapper<Cell>> cells)
{
	for (auto cell : cells) {
		neighbors_.push_back(cell);
	}
}

unsigned long Cell::get_x()
{
	return x_;
}

unsigned long Cell::get_y()
{
	return y_;
}
