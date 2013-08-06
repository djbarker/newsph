#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <list>
#include <boost/mpi/datatype.hpp>
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

	// more than checking equality this checks whether they are the same object
	bool is(Particle<Dim,TStep,NCol>&);

    //friend class boost::serialization::access;
	template<class Archive> void serialize(Archive& a, const unsigned int version);

	// properties
	size_t fluid;
	size_t wall;
	ParticleType type;
	nvect<Dim,quantity<position>>		pos[TStep];
	nvect<Dim,quantity<velocity>>		vel[TStep];
	nvect<Dim,quantity<acceleration>>	acc;
	quantity<IntDim<0,-(int)Dim,0>>			sigma;
	nvect<Dim,quantity<IntDim<0,-1,0>>>	gradC[NCol];

	typename std::list<Particle<Dim,TStep,NCol>*>::iterator lcg_position;
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

template<size_t Dim, size_t TStep, size_t NCol>
bool Particle<Dim,TStep,NCol>::is(Particle<Dim,TStep,NCol>& part)
{
	return this==&part;
}

} /* namespace sim */

// optimize sending via boost::mpi
namespace boost { namespace mpi {
  template <size_t Dim, size_t TStep, size_t NCol>
  struct is_mpi_datatype<sim::Particle<Dim,TStep,NCol> > : mpl::true_ { };
} }

#endif /* PARTICLE_H_ */
