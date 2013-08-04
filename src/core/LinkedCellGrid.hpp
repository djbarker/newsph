#ifndef LINKEDCELLGRID_HPP_
#define LINKEDCELLGRID_HPP_

#include <vector>
#include <list>
#include "Particle.hpp"
#include "../utils/utils.hpp"

namespace sim
{

enum PaddingLocation
{
	Top    = 0x1,  // 1
	Bottom = 0x2,  // 2
	Left   = 0x4,  // 4
	Right  = 0x8,  // 8
	Front  = 0x10, // 16
	Back   = 0x20, // 32
};

template<size_t Dim, typename T, size_t Padding=1>
class LinkedCellGrid
{
public:
	LinkedCellGrid(){};
	virtual ~LinkedCellGrid(){};

	void init(quantity<length> cell_size);

	Subscript<Dim> idxToSub(size_t idx);
	size_t idxToSub(Subscript<Dim> sub);
	std::list<T>& getCell(size_t idx);
	void place(std::list<T>& list);

	std::list<T> getPadding(size_t location);

private:
	quantity<length> cell_size;
	nvect<Dim,size_t> cell_counts;
	std::vector<std::list<T>> cells();
};

}

#endif /* LINKEDCELLGRID_HPP_ */
