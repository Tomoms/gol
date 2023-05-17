#include <algorithm>
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

std::vector<std::reference_wrapper<Cell>>& Cell::get_neighbors()
{
	return neighbors_;
}


unsigned long Cell::get_x() const
{
	return x_;
}

unsigned long Cell::get_y() const
{
	return y_;
}

bool Cell::is_alive() const
{
	return alive_;
}

bool Cell::is_alive_after_evolution() const
{
	unsigned char living_neighbors = std::count_if(neighbors_.begin(), neighbors_.end(), [](Cell& c) { return c.is_alive(); });
	return living_neighbors == 3 || (alive_ && living_neighbors == 2);
}

void Cell::set_alive(bool alive)
{
	alive_ = alive;
}

std::ostream& operator<<(std::ostream& os, Cell const & cell)
{
	os << "Cell x = " << cell.x_ << "; y = " << cell.y_ << "; is alive: " << cell.alive_ << "\n";
	os << "has got as neighbors the cells with coordinates:" << "\n";
	for (auto neighbor : cell.neighbors_) {
		os << "(" << neighbor.get().x_ << ", " << neighbor.get().y_ << ")" << "\n";
	}
	return os;
}
