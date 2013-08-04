#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>

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
	void serialize(std::ostream& out);
	void deserialize(std::istream& in);

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

template<size_t Dim, size_t TStep, size_t NCol>
void Particle<Dim,TStep,NCol>::serialize(std::ostream& out)
{
	out.write((char*)&fluid,sizeof(size_t));
	out.write((char*)&wall,sizeof(size_t));
	out.write((char*)&type,sizeof(ParticleType));

	for(size_t i=0; i<TStep; ++i)
	{
		pos[i].serialize(out);
		vel[i].serialize(out);
	}

	acc.serialize(out);
	out.write((char*)&sigma,sizeof(double));

	for(size_t i=0; i<NCol; ++i)
	{
		gradC[i].serialize(out);
	}
}

template<size_t Dim, size_t TStep, size_t NCol>
void Particle<Dim,TStep,NCol>::deserialize(std::istream& in)
{
	in.read((char*)&fluid,sizeof(size_t));
	in.read((char*)&wall,sizeof(size_t));
	in.read((char*)&type,sizeof(ParticleType));


	for(size_t i=0; i<TStep; ++i)
	{
		pos[i] = nvect<Dim,quantity<position>>().deserialize(in);
		vel[i] = nvect<Dim,quantity<velocity>>().deserialize(in);
	}

	acc = nvect<Dim,quantity<acceleration>>().deserialize(in);

	double sig_tmp;
	in.read((char*)&sig_tmp,sizeof(double));
	sigma = quantity<IntDim<0,-Dim,0>>(sig_tmp);

	for(size_t i=0; i<NCol; ++i)
	{
		gradC[i] = nvect<Dim,quantity<IntDim<0,-1,0>>>().deserialize(in);
	}
}

} /* namespace sim */
#endif /* PARTICLE_H_ */
