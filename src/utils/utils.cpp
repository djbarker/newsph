#include "utils.hpp"

using namespace std;

/*
 * Specialisations of template functions declared in utils.hpp
 */

namespace sim
{
namespace utils
{

template<>
size_t sub_to_idx<2>(const Subscript<2>& sub, const nvect<2,size_t>& extent)
{
	return sub[0]+sub[1]*((int)extent[0]);
}

template<>
size_t sub_to_idx<3>(const Subscript<3>& sub, const nvect<3,size_t>& extent)
{
	return sub[0]+((int)extent[0])*(sub[1]+sub[2]*((int)extent[1]));
}

template<>
Subscript<2> idx_to_sub<2>(size_t idx, const nvect<2,size_t>& extent)
{
	Subscript<2> out;
	out[0] = idx % extent[0];
	out[1] = (idx/extent[0]) % extent[1];
	return out;
}

template<>
Subscript<3> idx_to_sub<3>(size_t idx, const nvect<3,size_t>& extent)
{
	Subscript<3> out;
	out[0] = idx % extent[0];
	out[1] = (idx/extent[0]) % extent[1];
	out[2] = ((idx/extent[0])/extent[1]) % extent[2];
	return out;
}

template<>
nvect<2,size_t> calc_num_domains<2>(size_t nproc)
{
	size_t root = 1;
	for(size_t i=2;i<=sqrt(nproc);++i)
		if(nproc%i==0) root = i;

	return nvect<2,size_t>{root,nproc/root};
}

template<>
nvect<3,size_t> calc_num_domains<3>(size_t nproc)
{
 	size_t root1 = 1;
	for(size_t i=2;i<=pow(nproc,1./3.);++i)
		if(nproc%i==0) root1 = i;

	size_t root2 = 1;
	for(size_t i=2;i<=sqrt(nproc/root1);++i)
		if((nproc/root1)%i==0) root2 = i;

	return nvect<3,size_t>{root1,root2,nproc/(root1*root2)};
}

bool is_whitespace(char c)
{
	return c==' ' || c=='\t' || c=='\r' || c=='\n';
}

string strip_whitespace(string s)
{
	string out;
	for(auto c : s)
		if(!is_whitespace(c))
			out.push_back(c);
	return out;
}


} /* namespace utils */
} /* namespace sim */
