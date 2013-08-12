#include "Simulation.hpp"

#include <boost/mpi/nonblocking.hpp>
#include <boost/serialization/list.hpp>

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
				fluid_particles.push_back(part);
			}
		}

	cout << endl;
}

/*
 * Exchanges particles at the border with neighbouring processes and places
 * the received particles into the linked cell grid.
 */
template<>
void Simulation<2>::exchange()
{

	/*if(comm_rank==0)
	{
		for(int j=-1;j<(int)cells.cellCount()[1]+1;++j)
		{
			for(int i=-1;i<(int)cells.cellCount()[0]+1;++i)
			{
				bool out = false;
				auto tmp = Subscript<2>(i,j);
				if( _lcg_impl<2,particle_type,1,Bottom>::paddingMin(cells) <= tmp && tmp < _lcg_impl<2,particle_type,1,Bottom>::paddingMax(cells)  )
				{
					out = true;
					cout << "#";
				}

				if(!out)
				{
					cout << ".";
				}

				cout << "\t";
			}
			cout << endl;
		}
	}
	exit(0);*/

	/*
	 * For shifting the periods
	 */
	auto shift_period_up = [&](plist_type& l, size_t Dim)->void{
		for(auto& part : l)
		{
			part.pos[0][Dim] += gdomain.upper[Dim];
			part.pos[1][Dim] += gdomain.upper[Dim];
		}
	};

	auto shift_period_down = [&](plist_type& l, size_t Dim)->void{
		for(auto& part : l)
		{
			part.pos[0][Dim] -= gdomain.upper[Dim];
			part.pos[1][Dim] -= gdomain.upper[Dim];
		}
	};

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

	neighbour_particles.clear();

	/*
	 * Send & receive the data.
	 */

	// 8 sends and 8 receives
	mpi::request sreqs[8];
	mpi::request rreqs[8];
	plist_type sdata[8];
	plist_type rdata[8];

	Subscript<2> shifts[8] = { {-1, 0}, // left
							   { 1, 0}, // right
							   { 0, 1}, // top
							   { 0,-1}, // bottom
							   {-1,-1}, // bottom-left
							   {-1, 1}, // top-left
							   { 1, 1}, // top-right
							   { 1,-1}, // bottom-right
							 };

	// get the data from each part to send
	sdata[0] = cells.getBorder<Left>();
	sdata[1] = cells.getBorder<Right>();
	sdata[2] = cells.getBorder<Top>();
	sdata[3] = cells.getBorder<Bottom>();
	sdata[4] = cells.getBorder<Bottom|Left>();
	sdata[5] = cells.getBorder<Top|Left>();
	sdata[6] = cells.getBorder<Top|Right>();
	sdata[7] = cells.getBorder<Bottom|Right>();

	if(comm_rank==4)
		for(auto l : sdata)
			cout << l.size() << endl;

	size_t stags[8] = { Left, Right, Top, Bottom, Bottom|Left, Top|Left, Top|Right, Bottom|Right };
	size_t rtags[8] = { Right, Left, Bottom, Top, Top|Right, Bottom|Right, Bottom|Left, Top|Left }; // The process to our left (for example) is sending data to its right.

	// TODO: this could quite easily be generated to D dimensions as a member function for use in 3D
	auto do_swap = [&](size_t i)->void{
		// get the subscript & rank of the process we are sending to
		Subscript<2> dest_sub = utils::mod(domain_sub + shifts[i], domain_counts);
		int dest_rank = sub_to_idx(dest_sub,domain_counts);

		// if we are on a period boundary shift the data before we send it
		for(size_t d=0;d<2;++d)
		{
			if(dest_sub[d]!=domain_sub[d] && shifts[i][d]==-1 && domain_sub[d]==0)
				shift_period_up(sdata[i],d);

			if(dest_sub[d]!=domain_sub[d] && shifts[i][d]==+1 && domain_sub[d]==(int)domain_counts[d]-1)
				shift_period_down(sdata[i],d);
		}

		// send and receive
		sreqs[i] = comm.isend(dest_rank,stags[i],sdata[i]);
		rreqs[i] = comm.irecv(dest_rank,rtags[i],rdata[i]);
	};

	for(size_t i=0;i<8;++i)	do_swap(i);

	// wait for exchanges to finish
	mpi::wait_all(sreqs,sreqs+8);
	mpi::wait_all(rreqs,rreqs+8);

	/*
	 * Handle received data
	 */

	// move the received data into neighbour_particles
	for(plist_type& list : rdata)
	{
		neighbour_particles.splice(neighbour_particles.end(),list);
	}

	for(auto& p : neighbour_particles)
		p.type = GhostP; // DEBUG

	// add to linked cell grid
	cells.place(neighbour_particles,0);
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
