#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>
#include "../utils/utils.hpp"

namespace sim
{

enum ParticleType
{
	FluidP,
	WallP,
	GhostP,
	UnusedP
};

using namespace dims;

template<size_t Dim, size_t TStep, size_t NCol>
class Particle
{

public:
	Particle();
	virtual ~Particle(){};

    //friend class boost::serialization::access;
	template<class Archive> void serialize(Archive& a, const unsigned int version);

	// properties
	size_t fluid;
	size_t wall;
	ParticleType type;
	nvect<Dim,quantity<position>>		pos[TStep];
	nvect<Dim,quantity<velocity>>		vel[TStep];
	nvect<Dim,quantity<acceleration>>	acc;
	quantity<IntDim<0,-Dim,0>>			sigma;
	nvect<Dim,quantity<IntDim<0,-1,0>>>	gradC[NCol];
};

template<size_t Dim, size_t TStep, size_t NCol>
Particle<Dim,TStep,NCol>::Particle()
:fluid(0)
,wall(0)
,type(UnusedP)
{
}

template<size_t Dim, size_t TStep, size_t NCol> template<class Archive>
void Particle<Dim,TStep,NCol>::serialize(Archive& a, const unsigned int version)
{
	a & fluid;
	a & wall;
	a & type;
	a & pos; // boost automatically handles arrays
	a & vel;
	a & acc;
	a & sigma;
	a & gradC;
}

} /* namespace sim */
#endif /* PARTICLE_H_ */
