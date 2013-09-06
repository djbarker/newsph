#ifndef PARTICLEEXCEPTION_HPP_
#define PARTICLEEXCEPTION_HPP_

#include "../core/Particle.hpp"
#include <stdexcept>
#include <sstream>

template<class PType>
class ParticleException : public std::runtime_error
{
public:
	explicit ParticleException(const PType& part) noexcept;
	explicit ParticleException(const PType& part, std::string msg) noexcept;
	virtual ~ParticleException() noexcept;
	const char* what() const noexcept;


protected:
	const PType& _what;
	std::string msg;
};

template<class PType>
ParticleException<PType>::ParticleException(const PType& part) noexcept
:std::runtime_error("")
,_what(part)
,msg("")
{
}

template<class PType>
ParticleException<PType>::ParticleException(const PType& part, std::string msg) noexcept
:std::runtime_error("")
,_what(part)
,msg(msg)
{
}

template<class PType>
ParticleException<PType>::~ParticleException() noexcept
{
}

template<class PType>
const char* ParticleException<PType>::what() const noexcept
{
	using namespace std;

	stringstream sstr;
	sstr << msg << " due to " << _what;
	return sstr.str().c_str();
}

#endif /* PARTICLEEXCEPTION_HPP_ */
