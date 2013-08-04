#ifndef FLUID_H_
#define FLUID_H_

#include <iostream>
#include <dims.hpp>

namespace sim
{
	struct Fluid
	{
		bool gravity; // feels gravity?
		dims::quantity<dims::viscosity> viscosity; // dynamic viscosity
		dims::quantity<dims::density> density; // fluid density
		dims::quantity<dims::velocity> speed_of_sound;

		void serialize(std::ostream& out);
		void deserialize(std::istream& in);
	};

}

#endif /* FLUID_H_ */
