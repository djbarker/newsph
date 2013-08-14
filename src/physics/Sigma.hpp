#ifndef SIGMA_HPP_
#define SIGMA_HPP_

#include <dims.hpp>
#include "../kernels/ParticleDelta.hpp"

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
	void operator() (PType& a, PType& b, const kernels::ParticleDelta<Dim>& delta, Simulation<Dim>& sim)
	{
		//cout << W_ab << endl;
		a.sigma += delta.kernel;

		if(!a.is(b))
			b.sigma += delta.kernel;
	}

};

/*
 * Calculates the viscosity of the fluids particls
 */

template<int Dim, size_t Tstep>
struct ViscCalc {

	template<class PType>
	void operator() (PType& a, PType& b, const kernels::ParticleDelta<Dim>& delta, Simulation<Dim>& sim)
	{
		if(!a.is(b))
		{
			quantity<viscosity> mu_a = sim.fluidPhases()[a.fluid].viscosity;
			quantity<viscosity> mu_b = sim.fluidPhases()[a.fluid].viscosity;
			quantity<viscosity> mu = 2.0_number*mu_a*mu_b/(mu_a + mu_b);

			auto visc = mu*(pow<-2>(a.sigma) + pow<-2>(b.sigma))*delta.grad*(a.vel[Tstep]-b.vel[Tstep])/delta.dist;
			a.acc += visc/(sim.fluidPhases()[a.fluid].density*sim.parameters().V);
			b.acc -= visc/(sim.fluidPhases()[b.fluid].density*sim.parameters().V);
		}
	}
};

/*
 * Calculates the acceleration due to the pressure gradients
 */
template<int Dim>
struct GradPCalc {

	template<class PType>
	void operator() (PType& a, PType& b, const kernels::ParticleDelta<Dim>& delta, Simulation<Dim>& sim)
	{
		if(!a.is(b))
		{
			auto gradP = ((a.pressure*pow<-2>(a.sigma) + b.pressure*pow<-2>(b.sigma))*delta.grad*delta.unit);
			a.acc += gradP/(sim.fluidPhases()[a.fluid].density*sim.parameters().V);
			b.acc -= gradP/(sim.fluidPhases()[b.fluid].density*sim.parameters().V);
		}
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
		part.acc = make_vect<Dim,quantity<acceleration>>(0.0); // expands to (0.,0.,...) depending on Dim
	}
};

/*
 * Calculates a pressure based upon the Tait equation
 */
template<int Dim>
struct TaitEquation {

	Simulation<Dim>& sim;

	TaitEquation(Simulation<Dim>& s):sim(s){}

	template<class PType>
	void operator() (PType& part)
	{
		Fluid f = sim.fluidPhases()[part.fluid];
		part.pressure = (f.density*pow<2>(f.speed_of_sound)/7.0_number)*( pow<7>(part.density[0]/f.density) - 1.0_number );
	}

};

/*
 * Sets the fluid density from Sigma
 */
template<int Dim>
struct DensityCalc {

	Simulation<Dim>& sim;

	DensityCalc(Simulation<Dim>& s):sim(s){}

	template<class PType>
	void operator() (PType& part)
	{
		Fluid f = sim.fluidPhases()[part.fluid];
		part.density[0] = part.sigma*f.density*sim.parameters().V;
	}

};

}
}


#endif /* SIGMA_HPP_ */
