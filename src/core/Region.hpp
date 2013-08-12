#ifndef REGION_HPP_
#define REGION_HPP_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>
#include "../utils/utils.hpp"

namespace sim
{

using namespace utils;

template<size_t Dim>
struct Region
{
	Region();

	nvect<Dim,quantity<position>> upper;
	nvect<Dim,quantity<position>> lower;
	bool ellipse;

	bool inside(const nvect<Dim,quantity<length>>& position);

	template<class Archive> void serialize(Archive& a, const unsigned int version);
};

template<size_t Dim>
Region<Dim>::Region()
:ellipse(false)
{
	// lower & upper default to (0,0..)
}

template<size_t Dim>
bool Region<Dim>::inside(const nvect<Dim,quantity<length>>& pos)
{
	if(ellipse)
		return test_ellipse<Dim>(pos,upper,lower);
	else
		return (lower<=pos) && (pos<=upper);
}

template<size_t Dim> template<class Archive>
void Region<Dim>::serialize(Archive& a, const unsigned int version)
{
	a & upper;
	a & lower;
	a & ellipse;
}

} /* namespace sim */


#endif /* REGION_HPP_ */
