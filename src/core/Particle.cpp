#include "Particle.hpp"

std::ostream& operator<<(std::ostream& out, const sim::ParticleType& t)
{
	switch(t)
	{
	case sim::UnusedP:
		out << "UnusedP";
		break;
	case sim::FluidP:
		out << "FluidP";
		break;
	case sim::WallP:
		out << "WallP";
		break;
	case sim::GhostP:
		out << "GhostP";
		break;
	default:
		out << "UnknownP (" << t << ")";
		break;
	}
	return out;
}
