#ifndef WENDLANDQUINTIC_HPP_
#define WENDLANDQUINTIC_HPP_

#include "dims.hpp"

namespace sim
{
namespace kernels
{

using namespace dims;

template<int Dim>
struct WendlandQuintic
{
	static quantity<IntDim<0,-Dim,0>>   Kernel(quantity<length> r, quantity<length> h) {
		auto q = r/h;
		if(q<2.0_number)
		{
			auto x = (1.0_number-q*0.5_number);
			x *= x;
			x *= x;
			return C*x*(2.0_number*q+1.0_number)*pow<-(int)Dim>(h);
		}
		else return quantity<IntDim<0,-Dim,0>>(0.);
	}

	static quantity<IntDim<0,-Dim-1,0>> Grad(quantity<length> r, quantity<length> h) {
		auto q = r/h;
		if(q<2.0_number)
		{
			auto x = (1.0_number-q*0.5_number);
			x = x*x*x;
			return 5.0_number*C*pow<-Dim-1>(h)*q*x;
		}
		else return quantity<IntDim<0,-Dim-1,0>>(0.);
	}

private:
	static const quantity<number> C; // normalization constant
};

} /* namespace kernels */
} /* namespace sim */


#endif /* WENDLANDQUINTIC_HPP_ */
