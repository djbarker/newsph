#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <list>
#include <boost/mpi/datatype.hpp>
#include <boost/intrusive/list.hpp>
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
	Particle(const Particle<Dim,TStep,NCol>&);
	virtual ~Particle(){};

	Particle<Dim,TStep,NCol>&  operator= (const Particle<Dim,TStep,NCol>&);

	// more than checking equality this checks whether they are the same object
	bool is(Particle<Dim,TStep,NCol>&);

    //friend class boost::serialization::access;
	template<class Archive> void serialize(Archive& a, const unsigned int version);

	template<size_t D, size_t T, size_t C>
	friend std::ostream& operator<<(std::ostream& out, const Particle<D,T,C>& p);

	// properties
	size_t fluid;
	size_t wall;
	size_t id;
	ParticleType type;
	nvect<Dim,quantity<position>>		pos[TStep];
	nvect<Dim,quantity<velocity>>		vel[TStep];
	nvect<Dim,quantity<acceleration>>	acc;
	quantity<IntDim<0,-(int)Dim,0>>		sigma;
	quantity<dims::density>				density[TStep];
	quantity<dims::pressure>			pressure;
	nvect<Dim,quantity<IntDim<0,-1,0>>>	gradC[NCol];

	// hooks for multiply linked lists
	boost::intrusive::list_member_hook<>
		main_hook,			// for iterating over fluid/wall particles
		lcg_hook,			// for iterating over linked-cell-grid cells
		neighbour_hook;     // for iterating over particles either sent or received from neighbouring processes
};

template<class T> using MainList      = boost::intrusive::list<T, boost::intrusive::member_hook<T, boost::intrusive::list_member_hook<>, &T::main_hook> >;
template<class T> using LCGList       = boost::intrusive::list<T, boost::intrusive::member_hook<T, boost::intrusive::list_member_hook<>, &T::lcg_hook> >;
template<class T> using NeighbourList = boost::intrusive::list<T, boost::intrusive::member_hook<T, boost::intrusive::list_member_hook<>, &T::neighbour_hook> >;

template<size_t Dim, size_t TStep, size_t NCol>
Particle<Dim,TStep,NCol>::Particle()
:fluid(0)
,wall(0)
,id(std::numeric_limits<size_t>::max()-1)
,type(UnusedP)
{
}

template<size_t Dim, size_t TStep, size_t NCol>
Particle<Dim,TStep,NCol>::Particle(const Particle<Dim,TStep,NCol>& part)
:fluid(part.fluid)
,wall(part.wall)
,id(part.id)
,type(part.type)
,sigma(part.sigma)
,pressure(part.pressure)
{
	for(size_t t=0;t<TStep;++t)
	{
		pos[t] = part.pos[t];
		vel[t] = part.vel[t];
		density[t] = part.density[t];
	}

	for(size_t c=0;c<NCol;++c)
	{
		gradC[c] = part.gradC[c];
	}
}

template<size_t Dim, size_t TStep, size_t NCol>
Particle<Dim,TStep,NCol>& Particle<Dim,TStep,NCol>::operator=(const Particle<Dim,TStep,NCol>& part)
{
	if(&part!=this)
	{
		fluid(part.fluid);
		wall(part.wall);
		id(part.id);
		type(part.type);
		sigma(part.sigma);
		pressure(part.pressure);

		for(size_t t=0;t<TStep;++t)
		{
			pos[t] = part.pos[t];
			vel[t] = part.vel[t];
			density[t] = part.density[t];
		}

		for(size_t c=0;c<NCol;++c)
		{
			gradC[c] = part.gradC[c];
		}
	}
	return *this;
}

template<size_t Dim, size_t TStep, size_t NCol> template<class Archive>
void Particle<Dim,TStep,NCol>::serialize(Archive& a, const unsigned int version)
{
	a & fluid;
	a & wall;
	a & id;
	a & type;
	a & pos; // boost automatically handles arrays
	a & vel;
	a & acc;
	a & sigma;
	a & density;
	a & pressure;
	a & gradC;
}

template<size_t Dim, size_t TStep, size_t NCol>
bool Particle<Dim,TStep,NCol>::is(Particle<Dim,TStep,NCol>& part)
{
	return this==&part;
}

template<size_t Dim, size_t TStep, size_t NCol>
std::ostream& operator<<(std::ostream& out, const sim::Particle<Dim,TStep,NCol>& part)
{
	out << "sim::Particle<" << Dim << ',' << TStep << ',' << NCol << ">::{";
	out << " id=" << part.id << " | " << "type=" << part.type << " | ";
	for(size_t t=0;t<TStep;++t)
		out << "pos[" << t << "]=" << part.pos[t] << " | ";
	out << '}';
	return out;
}

} /* namespace sim */

// TODO: Currently we cannot optimize sending via boost::mpi since Particle is not a POD (also nvect is not POD)
/*namespace boost { namespace mpi {
  template <size_t Dim, size_t TStep, size_t NCol>
  struct is_mpi_datatype<sim::Particle<Dim,TStep,NCol> > : public mpl::true_ { };
} }*/


#endif /* PARTICLE_H_ */
