#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <vect.hpp>
#include <dims.hpp>
#include "../core/Particle.hpp"

using namespace dims;

/*
 * Convenience functions for working with indices/subscripts
 */

// note: subscripts can be negative, extents cannot
template<size_t Dim> using Subscript = nvect<Dim,int>;
template<size_t Dim> using Extent = nvect<Dim,size_t>;

template<size_t Dim> size_t	sub_to_idx(Subscript<Dim>& sub, Extent<Dim>& extent);

template<size_t Dim> Subscript<Dim> idx_to_sub(size_t idx, Extent<Dim>& extent);

/*
 * Decompose the domain
 */

template<size_t Dim> nvect<Dim,size_t> calc_num_domains(size_t nproc);

/*
 * Various other functions
 */

template<size_t Dim>
bool test_ellipse(const nvect<Dim,quantity<length>>& pos
		, const nvect<Dim,quantity<length>>& lower
		, const nvect<Dim,quantity<length>>& upper)
{
	auto axis = upper - lower;
	auto orig = (lower+upper)*quantity<number>(0.5);

	nvect<Dim,quantity<length>> test = (pos-orig)*(pos-orig)*quantity<number>(4.0)/axis;

	return test.sum() < quantity<length>(1.0);
}

template<size_t Dim, typename T2, typename T1>
nvect<Dim,T2> vector_to_nvect(const std::vector<T1>& v)
{
	nvect<Dim,T2> out;
	for(size_t i=0;i<Dim;++i)
		out[i] = T2(v[i]);
	return out;
}

bool is_whitespace(char c);
std::string strip_whitespace(std::string s);

// only works if a > -b
template<class T>
int mod(int a, T b)
{
	return (a>0?a%(int)b:a+(int)b);
}

#endif /* UTILS_HPP_ */
