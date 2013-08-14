#ifndef LINKEDCELLGRID_HPP_
#define LINKEDCELLGRID_HPP_

#include <vector>
#include <list>
#include <stdexcept>
#include <boost/pool/pool_alloc.hpp>
#include "Particle.hpp"
#include "Region.hpp"
#include "../utils/utils.hpp"

namespace sim
{

using namespace std;

enum PaddingLocation
{
	Top    = 0x1,  // 1
	Bottom = 0x2,  // 2
	Left   = 0x4,  // 4
	Right  = 0x8,  // 8
	Front  = 0x10, // 16
	Back   = 0x20, // 32
};

// forward declaration
template<size_t _Dim, class _T, int _Padding, size_t Loc> struct _lcg_impl;

template< size_t Dim, typename T, size_t Padding=1>
class LinkedCellGrid
{
public:

	// type which is returned when making copies of data
	typedef std::list<T,boost::fast_pool_allocator<T>> plist_type;

	LinkedCellGrid(){};
	virtual ~LinkedCellGrid(){};

	void init(qvect<Dim,length> cell_size, Extent<Dim> cell_counts, qvect<Dim,length> lower);

	Subscript<Dim> idxToSub(size_t idx);
	size_t subToIdx(const Subscript<Dim>& sub);
	Subscript<Dim> posToSub(const nvect<Dim,quantity<position>>&);
	std::list<T*>& getCell(size_t idx);
	Extent<Dim> cellCount() const;

	template<template<class U, class A> class Container, class A>
	void place(Container<T*,A>& list, size_t tstep);

	template<template<class U, class A> class Container, class A>
	void place(Container<T,A>& list, size_t tstep);

	void clear();

	// functions which are forwarded
	template<size_t Loc> plist_type	getBorder();
	template<size_t Loc> void		clearPadding();

private:
	template<size_t _Dim, class _T, int _Padding, size_t Loc> friend struct _lcg_impl;
	void copyCellContents(plist_type& out, const Subscript<Dim>& cell_sub);

	qvect<Dim,length> lower;
	qvect<Dim,length> cell_sizes;
	Extent<Dim> cell_counts; // including padding
	nvect<Dim,int> cell_counts_unpadded;
	std::vector<std::list<T*>> cells;
};

template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::init(qvect<Dim,length> cell_sizes, Extent<Dim> cell_counts, qvect<Dim,length> lower)
{
	this->cell_sizes = cell_sizes;
	this->cell_counts = cell_counts;
	this->cell_counts_unpadded = vect_cast<int>(cell_counts);
	this->lower = lower;

	// total number of cells
	size_t ncells = 1;
	for(size_t i=0;i<Dim;++i)
	{
		this->cell_counts[i] += 2*Padding;
		ncells *= this->cell_counts[i];
	}

	// create empty cells
	cells.resize(ncells);
}

/*
 * The convert the index of a cell to a subscript. The subscript does not include
 * padding so in 2D - for example - the first cell has subscript (-Paddding,-Padding)
 * but index 0.
 */
template<size_t Dim, typename T, size_t Padding>
Subscript<Dim> LinkedCellGrid<Dim,T,Padding>::idxToSub(size_t idx)
{
	return idx_to_sub(idx,cell_counts)-make_vect<Dim,int>(Padding);
}

/*
 * Convert a subscript to a cell index. This automatically accounts for padding.
 * So in 2D - for example - the subscript (-P,-P) returns the index zero, where
 * P is the number of padding cells.
 */
template<size_t Dim, typename T, size_t Padding>
size_t LinkedCellGrid<Dim,T,Padding>::subToIdx(const Subscript<Dim>& sub)
{
	return sub_to_idx(sub+make_vect<Dim,int>(Padding),cell_counts);
}

/*
 * Returns the correct subscript for a position. Note, if the position is
 * inside the lower-left padded region - say - then it will return (-1,-1).
 */
template<size_t Dim, typename T, size_t Padding>
Subscript<Dim> LinkedCellGrid<Dim,T,Padding>::posToSub(const nvect<Dim,quantity<position>>& pos)
{
	return vect_cast<int>(utils::discard_dims((pos-lower)/cell_sizes)); // "cast" to int
}

/*
 * Returns the list of particle pointers which corresponds to the given index.
 */
template<size_t Dim, typename T, size_t Padding>
std::list<T*>& LinkedCellGrid<Dim,T,Padding>::getCell(size_t idx)
{
	return cells[idx];
}

/*
 * Returns the unpadded cell count in each direction.
 */
template<size_t Dim, typename T, size_t Padding>
Extent<Dim> LinkedCellGrid<Dim,T,Padding>::cellCount() const
{
	return Extent<Dim>(cell_counts_unpadded);
}

// place particle_type
template<size_t Dim, typename T, size_t Padding>
template<template<class U, class A> class Container, class A>
void LinkedCellGrid<Dim,T,Padding>::place(Container<T,A>& parts, size_t tstep)
{
	for(T& t : parts)
	{
		/*if(subToIdx(posToSub(t.pos[tstep]))<0 || subToIdx(posToSub(t.pos[tstep]))>=cells.size())
		{
			cout << lower << ": invalid index " << subToIdx(posToSub(t.pos[tstep])) << " caused by sub " << posToSub(t.pos[tstep]) << " caused by position " << t.pos[tstep] << endl;
		}*/
		cells[subToIdx(posToSub(t.pos[tstep]))].push_back(&t);
	}
}

// place particle_type*
template<size_t Dim, typename T, size_t Padding>
template<template<class U, class A> class Container, class A>
void LinkedCellGrid<Dim,T,Padding>::place(Container<T*,A>& parts, size_t tstep)
{
	for(T* ptr : parts)
	{
		cells[subToIdx(posToSub(ptr->pos[tstep]))].push_back(ptr);
	}
}

/*
 * Appends the contents of the cell specified by a subscript to the list
 * given as the first paramter. Note that this copies the data.
 */
template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::copyCellContents(plist_type& out, const Subscript<Dim>& sub)
{
	size_t idx = subToIdx(sub);

	// copy data from T* pointers into plist_type object
	for(T* pT : cells[idx])
	{
		out.push_back(*pT);
	}
}

/*
 * Clears all particles from the linked cell grid.
 */
template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::clear()
{
	for(auto& cell : cells)
		cell.clear();
}



/*
 * Implementation Note:
 *
 * Really we would want to partially these member function for Dim and Loc
 * but we cannot partially specialize members without partially specializing
 * the whole class so we are forced to forward the implementation. To
 * an intermediate struct which we can partially specialize.
 *
 * This leads to more verbose syntax here but maintains compile time safety,
 * without affecting syntax elsewhere. E.g. if we accidently tried
 * to call getBorder<Left|Right>() it would not compile because
 * there is no specilization of _lcg_impl for Left|Right. (Though the compiler
 * message would be like trying to read Linear A)
 */

/*
 * Returns a copy of the data that is in the specified border region of the
 * linked cell grid. Returns a copy so we can use it as an boost::mpi buffer.
 */
template<size_t Dim, class T, size_t Padding> template<size_t Loc>
typename LinkedCellGrid<Dim,T,Padding>::plist_type LinkedCellGrid<Dim,T,Padding>::getBorder()
{
	plist_type out;

	// get limits of the border region cell subscripts
	Subscript<Dim> bmin = _lcg_impl<Dim,T,Padding,Loc>::borderMin(*this);
	Subscript<Dim> bmax = _lcg_impl<Dim,T,Padding,Loc>::borderMax(*this);

	// copy the cell contents
	utils::multi_for(bmin,bmax,[&](const Subscript<Dim>& loop_pos)->void{
		copyCellContents(out,loop_pos);
	});

	return out;
}

/*
 * Clears the padding at the specified location. Note that this doesn't delete
 * the particles it just removes them from the linked cell grid.
 */
template<size_t Dim, class T, size_t Padding> template<size_t Loc>
void LinkedCellGrid<Dim,T,Padding>::clearPadding()
{
	// get limits of the border region cell subscripts
	Subscript<Dim> min = _lcg_impl<Dim,T,Padding,Loc>::paddingMin(*this);
	Subscript<Dim> max = _lcg_impl<Dim,T,Padding,Loc>::paddingMax(*this);

	// copy the cell contents
	utils::multi_for(min,max,[&](const Subscript<Dim>& loop_pos)->void{
		cells[subToIdx(loop_pos)].clear();
	});

}

///////////////////////////////////////////////////////////////////////////////
//                                     2D                                    //
///////////////////////////////////////////////////////////////////////////////

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Left> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, 0 );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( Padding, lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( -Padding, 0 );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, lcg.cell_counts_unpadded[1] );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Right> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]-Padding, 0 );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], 0 );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]+Padding, lcg.cell_counts_unpadded[1] );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Bottom> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, 0 );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], Padding );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, -Padding );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], 0 );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Top> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, lcg.cell_counts_unpadded[1]-Padding );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], lcg.cell_counts_unpadded[1]+Padding );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Bottom|Left> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, 0 );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( Padding, Padding );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( -Padding, -Padding );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, 0 );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Top|Left> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, lcg.cell_counts_unpadded[1]-Padding );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( Padding, lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( -Padding, lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( 0, lcg.cell_counts_unpadded[1]+Padding );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Top|Right> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]-Padding, lcg.cell_counts_unpadded[1]-Padding );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], lcg.cell_counts_unpadded[1] );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]+Padding, lcg.cell_counts_unpadded[1]+Padding );
	}
};

template<class _T, int Padding>
struct _lcg_impl<2,_T,Padding,Bottom|Right> {
	static Subscript<2> borderMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]-Padding, 0 );
	}

	static Subscript<2> borderMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], Padding );
	}

	static Subscript<2> paddingMin(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0], -Padding );
	}

	static Subscript<2> paddingMax(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		return Subscript<2>( lcg.cell_counts_unpadded[0]+Padding, 0 );
	}
};


///////////////////////////////////////////////////////////////////////////////
//                                     3D                                    //
///////////////////////////////////////////////////////////////////////////////

}

#endif /* LINKEDCELLGRID_HPP_ */
