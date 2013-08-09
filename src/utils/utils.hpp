#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>
//#include <boost/serialization/string.hpp>
//#include <boost/serialization/version.hpp>
//#include <boost/serialization/utility.hpp>
#include <vect.hpp>
#include <dims.hpp>
#include "../core/Particle.hpp"

namespace sim
{
namespace utils
{

using namespace dims;

template<size_t N, class Dims> using qvect = nvect<N,dims::quantity<Dims>>;

/*
 * Convenience functions for working with indices/subscripts
 */

// note: subscripts can be negative, extents cannot
template<size_t Dim> using Subscript = nvect<Dim,int>;
template<size_t Dim> using Extent = nvect<Dim,size_t>;

// TODO: make these (and nvect) pass by rvalue
template<size_t Dim> size_t	sub_to_idx(Subscript<Dim> sub, Extent<Dim> extent);
template<size_t Dim> Subscript<Dim> idx_to_sub(size_t idx, Extent<Dim> extent);

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

// Vector version of discard_dims, which converts a nvect<N,quantity<Dim,T>> into an nvect<N,T>
template<size_t N, typename T, typename Dim>
nvect<N,T> discard_dims(const nvect<N,quantity<Dim,T>>& v)
{
	nvect<N,T> out;
	for(size_t i=0;i<N;++i)
		out[i] = discard_dims(v[i]);
	return out;
}


// Casts a general nvect<N,T> into one storing another type U.
// Note: T must be castable to U.
template<size_t N, typename T, typename U>
nvect<N,U> nvect_cast(const nvect<N,T>& v)
{
	nvect<N,U> out;
	for(size_t i=0;i<N;++i)
		out[i] = (U)v[i];
	return out;
}

// Casts an nvect<N,quantity<Dim,T>> to one storing "raw" data of type U
// i.e. nvect<N,U>
// Note: T must be castable to U.
template<size_t N, typename Dim, typename T, typename U>
nvect<N,U> qvect_cast(const nvect<N,quantity<Dim,T>>& v)
{
	nvect<N,U> out;
	for(size_t i=0;i<N;++i)
		out[i] = (U)discard_dims(v[i]);
	return out;
}

bool is_whitespace(char c);
std::string strip_whitespace(std::string s);

// only works if a > -b
template<class T, class U>
T mod(T a, U b)
{
	return (a>=0?a%(int)b:a+(int)b);
}

// vectorized version
template<size_t N, typename T, typename U>
nvect<N,T> mod(const nvect<N,T>& a, const nvect<N,U>& b)
{
	nvect<N,T> out;
	for(size_t i=0;i<N;++i)
		out[i] = mod(a[i],b[i]);
	return out;
}


} /* namespace utils */
} /* namespace sim */

/*
 * Functions for serializing quantities and nvects with boost::serialization.
 */

namespace boost { namespace serialization {

using namespace dims;

// load and save quantities
template<typename Dims, typename T, typename Archive>
void save(Archive& ar, const quantity<Dims,T>& qty, const unsigned int version)
{
	T t = discard_dims(qty);
	ar & t;
}

template<typename Dims, typename T, typename Archive>
void load(Archive& ar, quantity<Dims,T>& qty, const unsigned int version)
{
	T t;
	ar & t;
	qty = quantity<Dims,T>(t);
}

// since quantity is a class template we must type this out manually
// instead of using BOOST_SERIALIZATION_SPLIT_FREE
template<typename Dims, typename T, typename Archive>
inline void serialize(Archive & ar, quantity<Dims,T>& qty, const unsigned int file_version)
{
    split_free(ar, qty, file_version);
}

}}

// serialize nvects
template<size_t N, typename T, typename Archive>
void serialize(Archive& ar, nvect<N,T>& vect, const unsigned int version)
{
	for(size_t i=0;i<N;++i)
		ar & vect[i];
}

#endif /* UTILS_HPP_ */
