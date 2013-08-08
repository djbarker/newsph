#ifndef LINKEDCELLGRID_HPP_
#define LINKEDCELLGRID_HPP_

#include <vector>
#include <list>
#include <stdexcept>
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
template<size_t _Dim, class _T, size_t _Padding, size_t Loc> struct _lcg_impl;

template<size_t Dim, typename T, size_t Padding=1>
class LinkedCellGrid
{
public:
	LinkedCellGrid(){};
	virtual ~LinkedCellGrid(){};

	void init(qvect<Dim,length> cell_size, Extent<Dim> cell_counts, qvect<Dim,length> lower);

	Subscript<Dim> idxToSub(size_t idx);
	size_t subToIdx(Subscript<Dim> sub);
	Subscript<Dim> posToSub(const nvect<Dim,quantity<position>>&);
	std::list<T>& getCell(size_t idx);

	template<template<class U, class A> class Container, class U, class A>
	void place(Container<U,A>& list, size_t tstep);

	template<template<class U, class A> class Container, class A>
	void place(Container<T,A>& list, size_t tstep);

	void clear();

	// functions which are forwarded to _lcg_imp
	template<size_t Loc> std::list<T> getBorder();

private:

	template<size_t _Dim, class _T, size_t _Padding, size_t Loc> friend struct _lcg_impl;
	void appendCellContents(std::list<T>& out,Subscript<Dim> cell_sub);

	qvect<Dim,length> lower;
	qvect<Dim,length> cell_sizes;
	nvect<Dim,size_t> cell_counts;
	std::vector<std::list<T>> cells;
};

template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::init(qvect<Dim,length> cell_sizes, Extent<Dim> cell_counts, qvect<Dim,length> lower)
{
	this->cell_sizes = cell_sizes;
	this->cell_counts = cell_counts;
	this->lower = lower;

	// total number of cells
	size_t ncells = 1;
	for(size_t i=0;i<Dim;++i)
	{
		cell_counts[i] += 2*Padding;
		ncells *= cell_counts[i];
	}

	cout << ncells << endl;

	// create empty cells
	cells.resize(ncells);
}

template<size_t Dim, typename T, size_t Padding>
Subscript<Dim> LinkedCellGrid<Dim,T,Padding>::idxToSub(size_t idx)
{
	return idx_to_sub(idx,cell_counts)-Subscript<Dim>(Padding);
}

template<size_t Dim, typename T, size_t Padding>
size_t LinkedCellGrid<Dim,T,Padding>::subToIdx(Subscript<Dim> sub)
{
	return sub_to_idx(sub+Subscript<Dim>(Padding),cell_counts);
}

template<size_t Dim, typename T, size_t Padding>
Subscript<Dim> LinkedCellGrid<Dim,T,Padding>::posToSub(const nvect<Dim,quantity<position>>& pos)
{
	return Subscript<Dim>(discard_dims((pos-lower)/cell_sizes)); // "cast" to int
}

template<size_t Dim, typename T, size_t Padding>
std::list<T>& LinkedCellGrid<Dim,T,Padding>::getCell(size_t idx)
{
	return cells[idx];
}

// place particle_type
template<size_t Dim, typename T, size_t Padding>
template<template<class U, class A> class Container, class A>
void LinkedCellGrid<Dim,T,Padding>::place(Container<T,A>& parts, size_t tstep)
{
	for(T t : parts)
	{
		cells[subToIdx(posToSub(t->pos[tstep]))].push_back(t);
	}
}

// place particle_type*
template<size_t Dim, typename T, size_t Padding>
template<template<class U, class A> class Container, class U, class A>
void LinkedCellGrid<Dim,T,Padding>::place(Container<U,A>& parts, size_t tstep)
{
	for(U t : parts)
	{
		cells[subToIdx(posToSub(t.pos[tstep]))].push_back(&t);
	}
}

template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::appendCellContents(std::list<T>& out, Subscript<Dim> sub)
{
	size_t idx = subToIdx(sub);

	// copy cell[idx] contents to the end of out
	std::copy(cells[idx].begin(), cells[idx].end(),
			  std::back_insert_iterator<std::list<T> >(out));
}

// clears references to
template<size_t Dim, typename T, size_t Padding>
void LinkedCellGrid<Dim,T,Padding>::clear()
{
	// TODO: clear
}

/*
 * Really we would want to partially this member function for Dim and Loc
 * but we cannot partially specialize members without partially specializing
 * the whole class so we are forced to forward the implementation. To
 * an intermediate struct which we can partially specialize.
 *
 * This leads to more verbose syntax here but maintains compile time safety,
 * without affecting syntax elsewhere. E.g. if we accidently tried
 * to call getBorder<Left|Right>() it would not compile because
 * there is no specilization of _lcg_impl for Left|Right.
 */

template<size_t Dim, class T, size_t Padding> template<size_t Loc>
std::list<T> LinkedCellGrid<Dim,T,Padding>::getBorder()
{
	return _lcg_impl<Dim,T,Padding,Loc>::getBorder(*this);
}

// partial specializations for Dim and Loc

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Left> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int j=0;j<(int)lcg.cell_counts[1]-2*(int)Padding;++j)
			for(int i=0;i<(int)Padding;++i)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Right> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int j=0;j<(int)lcg.cell_counts[1]-2*(int)Padding;++j)
			for(int i=lcg.cell_counts[0]-3*Padding;i<(int)lcg.cell_counts[0]-2*(int)Padding;++i)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Top> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=0;i<(int)lcg.cell_counts[0]-2*(int)Padding;++i)
			for(int j=lcg.cell_counts[1]-3*Padding;j<(int)lcg.cell_counts[1]-2*(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Bottom> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=0;i<(int)lcg.cell_counts[0]-2*(int)Padding;++i)
			for(int j=0;j<(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Bottom|Left> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=0;i<(int)Padding;++i)
			for(int j=0;j<(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Bottom|Right> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=lcg.cell_counts[0]-3*Padding;i<(int)lcg.cell_counts[0]-2*(int)Padding;++i)
			for(int j=0;j<(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Top|Right> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=lcg.cell_counts[0]-3*Padding;i<(int)lcg.cell_counts[0]-2*(int)Padding;++i)
			for(int j=lcg.cell_counts[1]-3*Padding;j<(int)lcg.cell_counts[1]-2*(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};

template<class _T, size_t Padding>
struct _lcg_impl<2,_T,Padding,Top|Left> {
	static std::list<_T> getBorder(LinkedCellGrid<2,_T,Padding>& lcg)
	{
		std::list<_T> out;
		for(int i=0;i<(int)Padding;++i)
			for(int j=lcg.cell_counts[1]-3*Padding;j<(int)lcg.cell_counts[1]-2*(int)Padding;++j)
				lcg.appendCellContents(out,Subscript<2>(i,j));
		return out;
	}
};


}

#endif /* LINKEDCELLGRID_HPP_ */
