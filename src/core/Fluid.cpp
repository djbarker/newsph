#include "Fluid.h"

using namespace std;

namespace sim
{

void Fluid::serialize(ostream& out)
{
	out.write((char*)&gravity,sizeof(bool));
	out.write((char*)&viscosity,sizeof(double));
	out.write((char*)&density,sizeof(double));
	out.write((char*)&speed_of_sound,sizeof(double));
}

void Fluid::deserialize(istream& in)
{
	in.read((char*)&gravity,sizeof(bool));
	in.read((char*)&viscosity,sizeof(double));
	in.read((char*)&density,sizeof(double));
	in.read((char*)&speed_of_sound,sizeof(double));
}

}
