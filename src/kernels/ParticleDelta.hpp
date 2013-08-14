#ifndef PARTICLEDELTA_HPP_
#define PARTICLEDELTA_HPP_

#include "../utils/utils.hpp"

namespace sim
{
namespace kernels
{

using namespace dims;

/*
 * A class which stores the required values when calculating SPH sums,
 * such as distance, W, gradW etc.
 */
template<int Dim>
struct ParticleDelta
{
	quantity<length>			dist;	// the distance between two particles
	qvect<Dim,number>			unit;	// the unit vector from a to b
	quantity<IntDim<0,-Dim,0>> 	kernel; // W(|r_a-r_b|,h)
	quantity<IntDim<0,-Dim-1,0>> grad;	// del W
};



}
}
#endif /* PARTICLEDELTA_HPP_ */
