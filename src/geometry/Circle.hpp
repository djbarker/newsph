#ifndef CIRCLE_HPP_
#define CIRCLE_HPP_

#include "Geometric.hpp"

namespace sim
{

template<size_t Dim>
class Circle : public Geometric<Dim>
{
public:
	quantity<position,dvect<Dim>> centre;
	quantity<length> radius;

	void reflect(const quantity<position,dvect<Dim>>& pos, quantity<position,dvect<Dim>>& result, quantity<number>& volume_factor);
	bool intercept(const quantity<position,dvect<Dim>>& a, const quantity<position,dvect<Dim>>& b);
};

template<size_t Dim>
void Circle<Dim>::reflect(const quantity<position,dvect<Dim>>& pos, quantity<position,dvect<Dim>>& result, quantity<number>& volume_factor)
{
	quantity<length,dvect<Dim>> dx = pos - centre;
	quantity<length> x_rad = dx->magnitude();
	quantity<length> dist = radius - x_rad;

	quantity<number,dvect<Dim>> norm = dx / x_rad;

	volume_factor = pow<Dim-1>((radius+dist)/(radius-dist));

	result = pos+2.*norm*dist;
}

}


#endif /* CIRCLE_HPP_ */
