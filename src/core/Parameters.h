#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <iostream>
#include <dims.hpp>
#include <vect.hpp>

namespace sim
{

template<size_t Dim>
struct Parameters
{
	quantity<length> dx; // spatial discretization length
	quantity<length> h;  // smoothing length
	nvect<Dim,quantity<acceleration>> gravity;
	quantity<pressure> bkg_pressure;
	quantity<dims::time> tmax;
	quantity<dims::time> tout;
	quantity<dims::time> dt;

	void serialize(std::ostream& out)
	{
		out.write((char*)&dx,sizeof(double));
		out.write((char*)&h,sizeof(double));
		gravity.serialize(out);
		out.write((char*)&bkg_pressure,sizeof(double));
		out.write((char*)&tmax,sizeof(double));
		out.write((char*)&tout,sizeof(double));
		out.write((char*)&dt,sizeof(double));
	}

	void deserialize(std::istream& in)
	{
		in.read((char*)&dx,sizeof(double));
		in.read((char*)&h,sizeof(double));
		gravity.deserialize(in);
		in.read((char*)&bkg_pressure,sizeof(double));
		in.read((char*)&tmax,sizeof(double));
		in.read((char*)&tout,sizeof(double));
		in.read((char*)&dt,sizeof(double));
	}
};

}


#endif /* PARAMETERS_H_ */
