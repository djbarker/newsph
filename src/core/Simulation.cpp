#include "Simulation.hpp"

#include <boost/mpi/nonblocking.hpp>
#include <boost/mpi/skeleton_and_content.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>

using namespace std;
using namespace boost;

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

	for(size_t i=0;i<=nx;++i)
		for(size_t j=0;j<=ny;++j)
		{
			particle_type part;

			part.fluid = fluid;
			part.type = FluidP;
			part.pos[1] = part.pos[0] = nvect<2,quantity<number>>(i,j)*params.dx + region.lower;

			if(ldomain.inside(part.pos[0]))
			{
				particle_store.push_back(part);
				fluid_particles.push_back(particle_store.back());
			}
		}
}

template<> Subscript<2> Simulation<2>::shifts[8] = { {-1, 0}, // left
		   { 1, 0}, // right
		   { 0, 1}, // top
		   { 0,-1}, // bottom
		   {-1,-1}, // bottom-left
		   {-1, 1}, // top-left
		   { 1, 1}, // top-right
		   { 1,-1}, // bottom-right
		 };

template<> size_t Simulation<2>::send_tags[8] = { Left, Right, Top, Bottom, Bottom|Left, Top|Left, Top|Right, Bottom|Right };
template<> size_t Simulation<2>::recv_tags[8] = { Right, Left, Bottom, Top, Top|Right, Bottom|Right, Bottom|Left, Top|Left }; // The process to our left (for example) is sending data to its right.

/*
 * Exchanges particles at the border with neighbouring processes and places
 * the received particles into the linked cell grid.
 */
template<>
void Simulation<2>::exchangeFull()
{
	/*
	 * Clear previously any exchanged particles
	 */

	cells.clearPadding<Left>();
	cells.clearPadding<Right>();
	cells.clearPadding<Bottom>();
	cells.clearPadding<Top>();
	cells.clearPadding<Bottom|Left>();
	cells.clearPadding<Bottom|Right>();
	cells.clearPadding<Top|Right>();
	cells.clearPadding<Top|Left>();

	for(size_t i=0;i<hc_elements(2);++i)
	{
		recv_particles[i].clear();
		send_particles[i].clear();
	}

	/*
	 * Send & receive the data.
	 */

	mpi::request reqs[hc_elements(2)*2]; // send and receive

	// get the data from each part to send
	send_particles[0] = cells.getBorder<Left>();
	send_particles[1] = cells.getBorder<Right>();
	send_particles[2] = cells.getBorder<Top>();
	send_particles[3] = cells.getBorder<Bottom>();
	send_particles[4] = cells.getBorder<Bottom|Left>();
	send_particles[5] = cells.getBorder<Top|Left>();
	send_particles[6] = cells.getBorder<Top|Right>();
	send_particles[7] = cells.getBorder<Bottom|Right>();

	for(size_t i=0;i<hc_elements(2);++i)
	{
		// send and receive
		reqs[i] 				= comm.isend(dest_ranks[i],send_tags[i],send_particles[i]);
		reqs[i+hc_elements(2)]	= comm.irecv(dest_ranks[i],recv_tags[i],recv_particles[i]);
	};

	// wait for exchanges to finish
	mpi::wait_all(reqs,reqs+2*hc_elements(2));

	/*
	 * Handle received data
	 */

	for(size_t i=0;i<hc_elements(2);++i)
	{
		// if we are on a period boundary, shift the data we received across the period
		for(size_t d=0;d<2;++d)
		{
			// received from +Period
			if(dest_periods[i][d]==PeriodDirec::Positive)
				for(auto& part : recv_particles[i])
				{
					part.first.pos[0][d] -= gdomain.upper[d];
					part.first.pos[1][d] -= gdomain.upper[d];
				}

			// received from zero
			if(dest_periods[i][d]==PeriodDirec::Negative)
				for(auto& part : recv_particles[i])
				{
					part.first.pos[0][d] += gdomain.upper[d];
					part.first.pos[1][d] += gdomain.upper[d];
				}
		}
	}

	// add to linked cell grid
	for(auto& l : recv_particles)
	{
		for(auto& p : l)
		{
			cells.place(p.first,0);
			p.first.type = GhostP; // FOR DEBUG
		}
	}
}

/*
 * Doesn't exchange the particles themselves or place them into the lcg
 * it just updates the values on previously exchanged particles.
 */
template<>
void Simulation<2>::exchangeData()
{
	/*
	 * When we copied the neighbouring particles for sending before we also stored
	 * a pointer to the originating particle. Use this to copy the values and then
	 * exhange them.
	 */

	mpi::content sc[hc_elements(2)];
	mpi::content rc[hc_elements(2)];
	for(size_t i=0;i<hc_elements(2);++i)
	{
		for(auto& list : send_particles)
			for(auto& ppair : list)
				ppair.first = *ppair.second;

		sc[i] = mpi::get_content(send_particles[i]);
		rc[i] = mpi::get_content(recv_particles[i]);
		comm.isend(dest_ranks[i],send_tags[i],sc[i]);
		comm.irecv(dest_ranks[i],recv_tags[i],rc[i]);
	}
}

template<>
vector<Subscript<2>> Simulation<2>::getStencil()
{
	vector<Subscript<2>> out = { {+1,+1}, {+0,+1}, {-1,+1}, {-1,+0}, {+0,+0} };
	return out;
}

template<>
vector<Subscript<3>> Simulation<3>::getStencil()
{
	vector<Subscript<3>> out = { {-1,-1,-1}, {-1,+0,-1}, {-1,+1,-1},
								 {+0,-1,-1}, {+0,+0,-1}, {+0,+1,-1},
								 {+1,-1,-1}, {+1,+0,-1}, {+1,+1,-1},
								 {-1,-1,+0}, {-1,+0,+0}, {-1,+1,+0},
								 {+0,+1,+0}, {+0,+0,+0} };
	return out;
}

} /* namespace sim */
