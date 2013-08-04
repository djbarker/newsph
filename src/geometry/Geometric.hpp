#ifndef GEOMETRIC_HPP_
#define GEOMETRIC_HPP_

#include <dims.hpp>
#include <vect.hpp>

namespace sim
{

/* Abstract base class for Geometry objects */
template<size_t Dim>
class Geometric
{
public:
	Geometric();
	virtual ~Geometric() throw();

	size_t id;
	void reflect(const quantity<position,dvect<Dim>>& pos, quantity<position,dvect<Dim>>& result, quantity<number>& volume_factor) = 0;
	bool intercept(const quantity<position,dvect<Dim>>& a, const quantity<position,dvect<Dim>>& b) = 0;
};

};


#endif /* GEOMETRIC_HPP_ */
