#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <spud>
#include <boost/pool/pool_alloc.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "LinkedCellGrid.hpp"
#include "Parameters.h"
#include "Fluid.h"
#include "Region.hpp"
#include "../utils/utils.hpp"

namespace sim
{

template<size_t Dim>
class Simulation
{
public:
	Simulation();
	virtual ~Simulation(){};

	// I/O
	void loadConfigXML(std::string fname);
	void loadWall(std::string fname);
	void writeOutput();

	template<class Archive> void serialize(Archive& a, const unsigned int version);


	void floodFill(const Region<Dim>&, const nvect<Dim,quantity<position>>&, size_t fluid);

	void exchange();

	typedef sim::Particle<Dim,2,2> particle_type;

private:

	void init();

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
	std::list<particle_type,boost::fast_pool_allocator<particle_type>> fluid_particles;
	std::list<particle_type,boost::fast_pool_allocator<particle_type>> wall_particles;

	LinkedCellGrid<Dim,particle_type*> cells;
};

template<size_t Dim>
Simulation<Dim>::Simulation()
{
	// init MPI variables
	int tmp;
	MPI_Comm_size(MPI_COMM_WORLD,&tmp);
	comm_size = tmp;
	MPI_Comm_rank(MPI_COMM_WORLD,&tmp);
	comm_rank = tmp;
}

template<size_t Dim>
void Simulation<Dim>::init()
{
	using namespace std;

	// domain arrangements
	domain_counts = calc_num_domains<Dim>(comm_size);
	if(!comm_rank) cout << "Domain decomposition: " << domain_counts << endl;

	auto tmp2 = nvect<Dim,quantity<number>>(domain_counts); // convert from a nvect of size_t to nvect of quantity<number,size_t> s
	nvect<Dim,quantity<length>> dom_sizes = (gdomain.upper - gdomain.lower)/tmp2;

	domain_sub = idx_to_sub<Dim>((size_t)comm_rank,domain_counts); // get our position amongst the domain_counts

	// calculate lower domain position
	ldomain.lower = nvect<Dim,quantity<number>>(domain_sub)*dom_sizes;
	ldomain.upper = ldomain.lower + dom_sizes;
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

	// setup local domain_counts
	init();

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
	}

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
void Simulation<Dim>::writeOutput()
{
	// TODO: use boost::iostreams to compress output
	using namespace std;
	stringstream fname_str;
	fname_str << root << "." << comm_rank << ".dat";

	ofstream fout(fname_str.str(),ios::binary);
	if(!fout.is_open())
	{
		cerr << "Error opening output file on proc " << comm_rank << endl;
	}

	boost::archive::binary_oarchive oarch(fout);
	oarch << *this;

	fout.close();
}

} /* namespace sim */

#endif /* SIMULATION_HPP_ */
