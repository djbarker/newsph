#ifndef FLUID_H_
#define FLUID_H_

#include <iostream>
#include <dims.hpp>

namespace sim
{
	struct Fluid
	{
		bool gravity; // feels gravity?
		dims::quantity<dims::viscosity>	viscosity; // dynamic viscosity
		dims::quantity<dims::density>	density; // fluid density
		dims::quantity<dims::velocity>	speed_of_sound;

		template<class Archive> void serialize(Archive& ar, const unsigned int version);
	};

	template<class Archive>
	void Fluid::serialize(Archive& ar, const unsigned int version)
	{
		ar & gravity;
		ar & viscosity;
		ar & density;
		ar & speed_of_sound;
	}

}

#endif /* FLUID_H_ */
