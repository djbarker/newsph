#include "Simulation.hpp"

using namespace std;

/*
 * Specializations for 2D/3D specific functionality
 */

namespace sim
{

template<>
void Simulation<2>::floodFill(const Region<2>& region, const nvect<2,quantity<position>>& xflood, size_t fluid)
{
	// TODO: actually flood fill!

	size_t nx = discard_dims((region.upper[0]-region.lower[0])/params.dx);
	size_t ny = discard_dims((region.upper[1]-region.lower[1])/params.dx);

	for(size_t i=0;i<nx;++i)
		for(size_t j=0;j<ny;++j)
		{
			particle_type part;

			part.fluid = 0;
			part.type = FluidP;
			part.pos[0] = nvect<2,quantity<number>>(i,j)*params.dx + region.lower;

			if(ldomain.inside(part.pos[0]))
			{
				fluid_particles.push_back(part);
			}
		}
}

template<>
void Simulation<2>::exchange()
{

	// left
	list<particle_type*> data = cells.getPadding(Left);

	mod(domain_sub[0]-1,domain_counts[0]);

}

} /* namespace sim */
