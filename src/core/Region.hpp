#ifndef REGION_HPP_
#define REGION_HPP_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>
#include "../utils/utils.hpp"

namespace sim
{

template<size_t Dim>
struct Region
{
	Region();

	nvect<Dim,quantity<position>> upper;
	nvect<Dim,quantity<position>> lower;
	bool ellipse;

	bool inside(const nvect<Dim,quantity<length>>& position);

	void serialize(std::ostream& out);
	void deserialize(std::istream& out);
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
		return (lower<pos) && (pos<upper);
}

template<size_t Dim>
void Region<Dim>::serialize(std::ostream& out)
{
	upper.serialize(out);
	lower.serialize(out);
	out.write((char*)&ellipse,sizeof(bool));
}

template<size_t Dim>
void Region<Dim>::deserialize(std::istream& in)
{
	upper.deserialize(in);
	lower.deserialize(in);
	in.read((char*)&ellipse,sizeof(bool));
}

} /* namespace sim */


#endif /* REGION_HPP_ */
