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
	mpi::request reqs[16]; // 8 sends and 8 receives

	/*
	 * Send & receive the data.
	 */

	cout << "proc " << comm_rank << " here A" << endl;

	// left
	plist_type sdata_left = cells.getBorder<Left>();
	cout << "proc " << comm_rank << " here A.1" << endl;
	plist_type rdata_left;
	Subscript<2> dest_sub = domain_sub;
	dest_sub[0] = utils::mod(domain_sub[0]-1,domain_counts[0]);
	int dest_rank = sub_to_idx(dest_sub,domain_counts);
	cout << "proc " << comm_rank << " here A.2: " << dest_sub << " -> " << dest_rank << ",\t" << sdata_left.size() << endl;
	reqs[0] = comm.isend(dest_rank,Left,sdata_left);
	cout << "proc " << comm_rank << " here A.3" << endl;
	reqs[8] = comm.irecv(dest_rank,Right,rdata_left);

	cout << "proc " << comm_rank << " here B" << endl;

	// right
	plist_type sdata_right = cells.getBorder<Right>();
	plist_type rdata_right;
	dest_sub = domain_sub;
	dest_sub[0] = mod(domain_sub[0]+1,domain_counts[0]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[1] = comm.isend(dest_rank,Right,sdata_right);
	reqs[9] = comm.irecv(dest_rank,Left,rdata_right);

	cout << "proc " << comm_rank << " here B" << endl;

	// top
	plist_type sdata_top = cells.getBorder<Top>();
	plist_type rdata_top;
	dest_sub = domain_sub;
	dest_sub[1] = mod(domain_sub[1]+1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[2]  = comm.isend(dest_rank,Top,sdata_top);
	reqs[10] = comm.irecv(dest_rank,Bottom,rdata_top);

	cout << "proc " << comm_rank << " here C" << endl;

	// bottom
	plist_type sdata_bottom = cells.getBorder<Bottom>();
	plist_type rdata_bottom;
	dest_sub = domain_sub;
	dest_sub[1] = mod(domain_sub[1]-1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[3]  = comm.isend(dest_rank,Bottom,sdata_bottom);
	reqs[11] = comm.irecv(dest_rank,Top,rdata_bottom);

	cout << "proc " << comm_rank << " here D" << endl;

	// top-left
	plist_type sdata_tl = cells.getBorder<Top|Left>();
	plist_type rdata_tl;
	dest_sub = domain_sub;
	dest_sub[0] = mod(domain_sub[0]-1,domain_counts[0]);
	dest_sub[1] = mod(domain_sub[1]+1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[4]  = comm.isend(dest_rank,Top|Left,sdata_tl);
	reqs[12] = comm.irecv(dest_rank,Bottom|Right,rdata_tl);

	cout << "proc " << comm_rank << " here E" << endl;

	// top-right
	plist_type sdata_tr = cells.getBorder<Top|Right>();
	plist_type rdata_tr;
	dest_sub = domain_sub;
	dest_sub[0] = mod(domain_sub[0]+1,domain_counts[0]);
	dest_sub[1] = mod(domain_sub[1]+1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[5]  = comm.isend(dest_rank,Top|Right,sdata_tr);
	reqs[13] = comm.irecv(dest_rank,Bottom|Left,rdata_tr);

	cout << "proc " << comm_rank << " here F" << endl;

	// bottom-right
	plist_type sdata_br = cells.getBorder<Bottom|Right>();
	plist_type rdata_br;
	dest_sub = domain_sub;
	dest_sub[0] = mod(domain_sub[0]+1,domain_counts[0]);
	dest_sub[1] = mod(domain_sub[1]-1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[6]  = comm.isend(dest_rank,Bottom|Right,sdata_br);
	reqs[14] = comm.irecv(dest_rank,Top|Left,rdata_br);

	cout << "proc " << comm_rank << " here G" << endl;

	// bottom-left
	plist_type sdata_bl = cells.getBorder<Bottom|Left>();
	plist_type rdata_bl;
	dest_sub = domain_sub;
	dest_sub[0] = mod(domain_sub[0]-1,domain_counts[0]);
	dest_sub[1] = mod(domain_sub[1]-1,domain_counts[1]);
	dest_rank = sub_to_idx(dest_sub,domain_counts);
	reqs[7]  = comm.isend(dest_rank,Bottom|Left,sdata_bl);
	reqs[15] = comm.irecv(dest_rank,Top|Right,rdata_bl);

	cout << "proc " << comm_rank << " here H" << endl;

	// wait for exchange to finish
	mpi::wait_all(reqs,reqs+16);

	cout << "proc " << comm_rank << " here I" << endl;

	/*
	 * Put particles into cells
	 */

	cells.place(rdata_top,0);
	cells.place(rdata_bottom,0);
	cells.place(rdata_left,0);
	cells.place(rdata_right,0);
	cells.place(rdata_tl,0);
	cells.place(rdata_tr,0);
	cells.place(rdata_br,0);
	cells.place(rdata_bl,0);

	cout << "proc " << comm_rank << " here J" << endl;
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
