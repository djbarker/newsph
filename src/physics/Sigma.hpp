#ifndef SIGMA_HPP_
#define SIGMA_HPP_

#include <dims.hpp>

namespace sim
{
namespace physics
{

using namespace dims;

/*
 * Adds particles' contributions to sigma
 */
template<size_t Dim>
struct SigmaCalc {

	template<class PType>
	void operator() (PType& a, PType& b, quantity<IntDim<0,-Dim,0>> W_ab, quantity<IntDim<0,-Dim-1,0>> gradW_ab)
	{
		a.sigma += W_ab;

		if(!a.is(b))
			b.sigma += W_ab;
	}

};

/*
 * Resets values at the start of an iteration to zero.
 */
template<size_t Dim>
struct ResetVals {

	template<class PType>
	void operator() (PType& part) {
		part.sigma(0.0);
		part.acc(0.0); // expands to (0.,0.,...) depending on Dim
	}
};


}
}


#endif /* SIGMA_HPP_ */
