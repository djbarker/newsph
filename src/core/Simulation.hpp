#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <initializer_list>
#include <spud>
#include <boost/pool/pool_alloc.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpi/communicator.hpp>
#include "LinkedCellGrid.hpp"
#include "Parameters.h"
#include "Fluid.h"
#include "Region.hpp"
#include "../utils/utils.hpp"

namespace sim
{

using namespace std;

template<size_t Dim>
class Simulation
{
public:
	typedef sim::Particle<Dim,2,2> particle_type;
	typedef std::list<particle_type,boost::fast_pool_allocator<particle_type>> plist_type;

	Simulation();
	virtual ~Simulation(){};

	// I/O
	void loadConfigXML(std::string fname);
	void loadWall(std::string fname);
	void writeOutput(size_t file_no);
	template<class Archive> void serialize(Archive& a, const unsigned int version);

	// Setup
	void init();
	void floodFill(const Region<Dim>&, const nvect<Dim,quantity<position>>&, size_t fluid);
	void assignParticleIds();

	// Simulate
	void exchange();
	void placeParticlesIntoLinkedCellGrid(size_t tstep);
	template<template<int> class K, typename... Fs> void doSPHSum(size_t tstep, Fs&&... fs);
	template<typename... Fs> void applyFunctions(Fs&&... fs);

	const Parameters<Dim>& parameters() const;
	const plist_type& fluidParticles() const;
	const plist_type& wallParticles() const;

private:

	std::vector<Subscript<Dim>> getStencil();

	boost::mpi::communicator comm;
	size_t comm_size;
	size_t comm_rank;

	Extent<Dim>    domain_counts; // how many domains in each dimension
	Subscript<Dim> domain_sub;    // local domain "subscript" within all procs
	Region<Dim>	   gdomain;       // global domain extent
	Region<Dim>	   ldomain;       // local domain extent

	std::string			root;   // output filename root
	Parameters<Dim>		params; // physical parameters
	std::vector<Fluid>	fluids; // fluid parameters

	// lists for particles
	plist_type fluid_particles;
	plist_type wall_particles;

	// TODO: perhaps make LCG hold iterators to the plist_type rather than pointers.
	LinkedCellGrid<Dim,particle_type*> cells;
};

template<size_t Dim>
Simulation<Dim>::Simulation()
{
	// init MPI variables
	comm_size = comm.size();
	comm_rank = comm.rank();
}

template<size_t Dim>
void Simulation<Dim>::init()
{
	using namespace std;

	// TODO: INIT LINKED CELL GRID AND ENSURE LOCAL DOMAINS ALIGN PROPERLY WITH EDGES

	// get cell sizes in each dimension (can be slightly off 2h to ensure they fit exactly in the domain)
	qvect<Dim,number> gnum_cells = gdomain.upper / (2.0_number*params.h);
	for(size_t i=0;i<Dim;++i) gnum_cells[i] = floor(gnum_cells[i]);
	qvect<Dim,length> cell_sizes = gdomain.upper / gnum_cells;
	if(!comm_rank) cout << "Cell sizes: " << cell_sizes/params.h << " * h" << endl;

	// domain arrangements
	domain_counts = calc_num_domains<Dim>(comm_size);
	if(!comm_rank) cout << "Domain decomposition: " << domain_counts << endl;

	Extent<Dim> global_cell_counts(discard_dims(gnum_cells)); // "cast" to size_t

	// calculate the number of cells for each processor in each dimension
	nvect<Dim,size_t> lnum_cells[comm_size];
	for(size_t p=0;p<comm_size;++p)
		for(size_t i=0;i<Dim;++i)
			if(p<floor(global_cell_counts[i]/domain_counts[i]))
				lnum_cells[p][i] = floor(global_cell_counts[i]/domain_counts[i]);
			else
				lnum_cells[p][i] = floor(global_cell_counts[i]/domain_counts[i])+1;

	// get the size in each dimension of our domain
	nvect<Dim,quantity<length>> dom_sizes = qvect<Dim,number>(lnum_cells[comm_rank])*cell_sizes;

	// get our position amongst the domain_counts
	domain_sub = idx_to_sub<Dim>((size_t)comm_rank,domain_counts);

	// calculate local domain physical position
	for(size_t i=0;i<Dim;++i)
		for(size_t j=0;j<(size_t)domain_sub[i];++j)
		{
			auto tmp = domain_sub;
			tmp[i] = j;
			ldomain.lower[i] += quantity<number>(lnum_cells[sub_to_idx<Dim>(tmp,domain_counts)][i])*cell_sizes[i];
		}

	ldomain.upper = ldomain.lower + dom_sizes;

	// initialize the linked cell grid.
	cells.init(cell_sizes,lnum_cells[comm_rank],ldomain.lower);

	cout << "P" << comm_rank <<" : " << ldomain.lower << "->" << ldomain.upper << endl;
}

template<size_t Dim>
void Simulation<Dim>::loadConfigXML(std::string fname)
{
	using namespace std;
	using namespace Spud;

	if(SPUD_NO_ERROR!=load_options(fname))
	{
		throw("Error loading XML options file!");
	}

	string tmps;
	vector<double> tmpvd;
	double tmpd;
	int tmpi;

	// check dimensions match
	int dim_opts;
	get_option("/geometry/dimension",dim_opts);

	if(dim_opts!=Dim)
	{
		if(!comm_rank) cerr << "Code not compiled for " << dim_opts << " dimensions." << endl
							<< "Please recompile with Dim=" << dim_opts << "." << endl;

		throw runtime_error("Incorrect dimensions!");
	}

	// get simulation domain
	vector<double> period;
	get_option("/geometry/period",period);
	gdomain.upper = vector_to_nvect<Dim,quantity<position>>(period);

	get_option("/file_io/root",root,"out");

	// load walls
	if(have_option("/file_io/walls"))
	{
		if(!comm_rank) cout << "Loading walls." << endl;
		get_option("/file_io/walls",tmps);
		loadWall(tmps);
	}

	// load fluids
	for(int i=0;i<option_count("/physics/fluid");++i)
	{
		stringstream sstr;
		sstr << "/physics/fluid[" << i << "]";
		string path = sstr.str();

		get_option(path+"/name",tmps);
		if(!comm_rank) cout << tmps << endl;

		Fluid tmpf;
		tmpf.gravity = have_option(path+"/gravity");
		get_option(path+"/density",tmpd);
		tmpf.density = quantity<density>(tmpd);

		if(have_option(path+"viscosity/dynamic"))
		{
			get_option(path+"viscosity/dyanmic",tmpd);
			tmpf.viscosity = quantity<viscosity>(tmpd);
		}
		else
		{
			get_option(path+"viscosity/kinematic",tmpd);
			tmpf.viscosity = quantity<IntDim<0,2,-1>>(tmpd)*tmpf.density; // convert to dynamic viscosity
		}
	}

	if(have_option("/physics/gravity"))
	{
		get_option("/physics/gravity",tmpvd);
		params.gravity = vector_to_nvect<Dim,quantity<acceleration>>(tmpvd);
	}

	get_option("/physics/atmospheric_pressure",tmpd,0.0);
	params.bkg_pressure = quantity<pressure>(tmpd);

	// SPH options
	double c0;
	get_option("/sph/c0",c0);

	int ref;
	get_option("sph/reference_fluid",ref,0);

	if(ref<0)
	{
		if(!comm_rank) cerr << "Reference fluid must be zero or above!";
		throw runtime_error("Invalid reference fluid!");
	}

	// set matching speeds of sound
	for(auto& fluid : fluids)
	{
		fluid.speed_of_sound = quantity<velocity>(c0)*sqrt(fluids[ref].density/fluid.density);
	}

	double hfac;
	get_option("/sph/h_factor",hfac);
	if(have_option("/sph/resolution/dx"))
	{
		get_option("/sph/resolution/dx",tmpd);
		params.dx = quantity<length>(tmpd);
		params.h = params.dx*quantity<number>(hfac);
	}
	else
	{
		get_option("/sph/resolution/h",tmpd);
		params.h = quantity<length>(tmpd);
		params.dx = params.h/quantity<number>(hfac);
	}

	double time;
	get_option("/time/t_max",time);
	params.tmax = quantity<dims::time>(time);
	get_option("/time/dt_write",time);
	params.tout = quantity<dims::time>(time);
	get_option("/time/dt_max",time);
	params.dt = quantity<dims::time>(time);

	// setup local domain_counts, linked cell grid, etc
	init();

	// perform any flood filling requested

	for(int i=0;i<option_count("/flood_fill");++i)
	{
		stringstream sstr;
		sstr << "/flood_fill[" << i << "]";
		string path = sstr.str();

		get_option(path+"/name",tmps);
		if(!comm_rank) cout << "Filling region \'" << tmps << "\'" << endl;

		Region<Dim> fill_region;

		get_option(path+"/fill_region/lower",tmpvd);
		fill_region.lower = vector_to_nvect<Dim,quantity<position>>(tmpvd);
		get_option(path+"/fill_region/upper",tmpvd);
		fill_region.upper = vector_to_nvect<Dim,quantity<position>>(tmpvd);
		fill_region.ellipse = have_option(path+"/fill_region/ellipse");

		// check region is valid
		if(fill_region.upper<fill_region.lower)
		{
			throw runtime_error("Flood fill region is invalid!");
		}

		get_option(path+"/start_point",tmpvd);
		auto start_point = vector_to_nvect<Dim,quantity<position>>(tmpvd);

		get_option(path+"/fluid",tmpi,0);

		floodFill(fill_region,start_point,(size_t)tmpi);

		if(!comm_rank) cout << "Done" << endl;
	}

	// finished setting up - assign ids to the particles
	assignParticleIds();

}

// TODO: should really use a proper parsing system (e.g boost::spirit)
template<size_t Dim>
void Simulation<Dim>::loadWall(std::string fname)
{
	using namespace std;

	/*
	 * Types used for parsing
	 */
	enum ObjectType {
		Circle,
		Line,
	};

	struct Object
	{
		ObjectType type;
		vector<pair<string,string>> properties;
	};

	/*
	 * Function to parse an object
	 */
	auto parse_object = [](ifstream& in) -> Object
	{
		Object object;

		// objects start with word "object"
		string token;
		getline(in,token,' ');
		if(token.compare("object")!=0)
			throw runtime_error("Unexpected keyword, expected \'object\'!");

		// object type
		getline(in,token,' ');

		     if(token.compare("circle")==0)	object.type = Circle;
		else if(token.compare("line")==0)	object.type = Line;
		else throw runtime_error("Unknown object type!");

		// "with" keyword
		getline(in,token,' ');
		if(token.compare("with")!=0)
			throw runtime_error("Expected \'with\' keyword after object type!");

		// now the properties
		string prop_name,prop_val;
		while(!in.eof())
		{
			// get name
			prop_name.clear();
			while(!in.eof() && prop_name.back()!='=') prop_name.push_back(in.get());

			// clear equals and whitespace
			prop_name.pop_back();
			prop_name = strip_whitespace(prop_name);

			// get value (as string)
			prop_val.clear();
			while(!in.eof() && (prop_val.back()!='&' && prop_val.back()!=';')) prop_val.push_back(in.get());
			char c = prop_val.back();

			// clear ampersand or semi-colon and  whitespace
			prop_val.pop_back();
			prop_val = strip_whitespace(prop_val);

			// store the property/value pair
			object.properties.push_back(make_pair(prop_name,prop_val));

			// semi-colon indicates end of object
			if(c==';') break;
		}

		if(is_whitespace(in.peek())) in.get();

		return object;
	};

	/*
	 * Open the file and parse
	 */
	ifstream fin(fname);
	if(!fin.is_open())
		throw runtime_error("Unable to open wall file!");

	while(!fin.eof())
	{
		try
		{
			Object object = parse_object(fin);

			switch(object.type)
			{
			case Circle:
				//cout << "Properties: " << endl;
				//for(auto p : object.properties)
				//	cout << p.first << ":" << p.second << endl;
				break;
			case Line:
				break;
			}
		}
		catch(runtime_error& e)
		{
			if(!comm_rank) cerr << "Unable to parse line in .wob file:" << endl
					 	 	    << "\t" << e.what() << endl;
		}
	}
}

template<size_t Dim>
void Simulation<Dim>::assignParticleIds()
{
	size_t gid = 0;

	if(comm_rank>0)
		comm.recv(comm_rank-1,comm_rank,gid);

	for(auto& part : fluid_particles)
	{
		part.id = gid;
		gid++;
	}

	for(auto& part : wall_particles)
	{
		part.id = gid;
		gid++;
	}

	if(comm_rank<comm_size-1)
		comm.send(comm_rank+1,comm_rank+1,gid);
}

template<size_t Dim>
const typename Simulation<Dim>::plist_type& Simulation<Dim>::fluidParticles() const
{
	return fluid_particles;
}

template<size_t Dim>
const typename Simulation<Dim>::plist_type& Simulation<Dim>::wallParticles() const
{
	return wall_particles;
}

template<size_t Dim> template<typename Archive>
void Simulation<Dim>::serialize(Archive& a, const unsigned int version)
{
	a & params;
	a & gdomain;
	a & ldomain;
	a & fluids;
	a & fluid_particles;
	a & wall_particles;
}

//BOOST_CLASS_VERSION(Simulation<2>,0)
//BOOST_CLASS_VERSION(Simulation<3>,0)

template<size_t Dim>
void Simulation<Dim>::writeOutput(size_t file_number)
{
	// TODO: use boost::iostreams to compress output
	using namespace std;
	stringstream fname_str;
	fname_str << root << '.' << file_number << '.' << comm_rank << ".dat";

	ofstream fout(fname_str.str(),ios::binary);
	if(!fout.is_open())
	{
		cerr << "Error opening output file on proc " << comm_rank << endl;
	}

	boost::archive::binary_oarchive oarch(fout);
	oarch << *this;

	fout.close();
}

template<size_t Dim>
const Parameters<Dim>& Simulation<Dim>::parameters() const
{
	return params;
}

/*
 * This function puts wall_particles and fluid_particles into the correct sublists based
 * upon their positions at the specified timestep.
 */
template<size_t Dim>
void Simulation<Dim>::placeParticlesIntoLinkedCellGrid(size_t tstep)
{
	cells.clear();

	cells.place(fluid_particles,tstep);
	cells.place(wall_particles,tstep);
}

/*
 * This function is used to actually perform the SPH sums over fluid & wall
 * particles. It accepts any callable objects of the form
 *
 * void func(particle_type& a, particle_type& b, quantity<IntDim<0,-Dim,0>> W_ab, quantity<IntDim<0,-Dim-1,0>> gradW_ab)
 *
 * Note; if a value is returned it is discarded.
 */
template<size_t Dim>
template<template<int> class Kernel, typename... Fs>
void Simulation<Dim>::doSPHSum(size_t tstep, Fs&&... fs)
{
	static_assert(sizeof...(Fs)>0,"No operations passed to doSPHSum()!");

	auto neighbour_cells = getStencil();

	cout << "HERE 0 doSPHSum" << endl;

	typename std::list<particle_type>::iterator itr = fluid_particles.begin();
	while(true)
	{
		Subscript<Dim> x_sub = cells.posToSub(itr->pos[tstep]);

		// iterate over nearby particles
		for(Subscript<Dim> dcell : neighbour_cells)
		{
			for(particle_type* part_b : cells.getCell(cells.subToIdx(x_sub+dcell)))
			{
				cout << part_b->pos[tstep] << endl;
				quantity<length>       dist_ab = (itr->pos[tstep]-part_b->pos[tstep]).magnitude();
				quantity<IntDim<0,-(int)Dim,0>>  W_ab = Kernel<2>::Kernel(dist_ab,params.h);
				quantity<IntDim<0,-1-(int)Dim,0>> dW_ab = Kernel<2>::Grad(dist_ab,params.h);

				cout << "HERE 2 doSPHSum" << endl;

				// for explanation of this line see: http://stackoverflow.com/questions/18077259/variadic-function-accepting-functors-callable-objects
				auto dummylist = { ((void)std::forward<Fs>(fs)(*itr,*part_b,W_ab,dW_ab),0)... };
				(void)dummylist; // stop the compiler warning about unused variable
			}

		}

		++itr;

		if(itr==fluid_particles.end())
			itr = wall_particles.begin();

		if(itr==wall_particles.end())
			break;
	}

}

/*
 * This funciton is used to apply functions / transformations to fluid and wall
 * particles. It accepts any callable objects of the form
 *
 * void func(particle_type& a)
 *
 * Note; if a value is returned it is discarded.
 */
template<size_t Dim>
template<typename... Fs>
void Simulation<Dim>::applyFunctions(Fs&&... fs)
{
	static_assert(sizeof...(Fs)>0,"No operations passed to applyFunctions()!");

	typename std::list<particle_type>::iterator itr = fluid_particles.begin();
	while(true)
	{
		// for explanation of this line see: http://stackoverflow.com/questions/18077259/variadic-function-accepting-functors-callable-objects
		auto dummylist = { ((void)std::forward<Fs>(fs)(*itr),0)... };
		(void) dummylist; // hide warning about unused variable

		++itr;

		if(itr==fluid_particles.end())
			itr = wall_particles.begin();

		if(itr==wall_particles.end())
			break;
	}
}

} /* namespace sim */

#endif /* SIMULATION_HPP_ */
