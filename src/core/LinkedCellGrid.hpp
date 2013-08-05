#ifndef LINKEDCELLGRID_HPP_
#define LINKEDCELLGRID_HPP_

#include <vector>
#include <list>
#include <stdexcept>
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

/*
 * Really we would want to partially specialize on Dim (and possibly Padding)
 * but we cannot partially specialize members without partially specializing
 * the whole class so we are forced to use ifs. This removes some compile-
 * time safety but it's more acceptable than the template work-arounds. It
 * would also be safer to have loc as a template parameter so we could tell
 * at compile time if we were accidently calling say getPadding<Front> with
 * Dims==2 which would be invalid.
 */

template<size_t Dim, class T, size_t Padding>
std::list<T> LinkedCellGrid<Dim,T,Padding>::getPadding(size_t loc)
{
	if(Dim==2)
	{
		// TODO: padding
		switch(loc)
		{
		case Top:
			break;
		case Bottom:
			break;
		case Left:
			break;
		case Right:
			break;
		case Top|Left:
			break;
		case Top|Right:
			break;
		case Bottom|Left:
			break;
		case Bottom|Right:
			break;
		default:
			throw std::logic_error("invald padding location for 2D LinkedCellGrid!"); // see above
		}
	}
	else if(Dim==3)
	{

	}

}

}

#endif /* LINKEDCELLGRID_HPP_ */
