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
template<int Dim>
struct SigmaCalc {

	template<class PType>
	void operator() (PType& a, PType& b, quantity<IntDim<0,-Dim,0>> W_ab, quantity<IntDim<0,-Dim-1,0>> gradW_ab)
	{
		cout << W_ab << endl;
		a.sigma += W_ab;

		if(!a.is(b))
			b.sigma += W_ab;
	}

};

/*
 * Resets values at the start of an iteration to zero.
 */
template<int Dim>
struct ResetVals {

	template<class PType>
	void operator() (PType& part)
	{
		part.sigma = quantity<IntDim<0,-Dim,0>>(0.0);
		part.acc = nvect<+Dim,quantity<acceleration>>(0.0); // expands to (0.,0.,...) depending on Dim
	}
};


}
}


#endif /* SIGMA_HPP_ */
