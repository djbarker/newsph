#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>
#include "../utils/utils.hpp"

namespace sim
{

template<size_t Dim>
struct Parameters
{
	quantity<length> dx; // spatial discretization length
	quantity<length> h;  // smoothing length
	nvect<Dim,quantity<acceleration>> gravity;
	quantity<pressure> bkg_pressure;
	quantity<dims::time> tmax;
	quantity<dims::time> tout;
	quantity<dims::time> dt;

	// convenience quantities
	quantity<IntDim<0,Dim,0>> V;

	template<typename Archive>
	void serialize(Archive& a, const unsigned int version)
	{
		a & dx;
		a & h;
		a & gravity;
		a & bkg_pressure;
		a & tmax;
		a & tout;
		a & dt;
		a & V;
	}
};

}


#endif /* PARAMETERS_H_ */
