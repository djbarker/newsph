#ifndef LISTSERIALIZER_HPP_
#define LISTSERIALIZER_HPP_

#include <boost/serialization/serialization.hpp>
#include <memory>

/*
 * A class which allows us to serialize boost::intrusive::list objects
 * by wrapping them in an object which knows about the memory allocation
 * patterns.
 */
template<template<class T_> class List, template<class T_, class Alloc_> class Store, class T, class Alloc=std::allocator<T> >
class ListSerializer
{
public:
	ListSerializer() = default;

	ListSerializer(Store<T,Alloc>& s, List<T>& l)
	:store(&s)
	,list(&l)
	{
	}

	template<class Archive>
	void save(Archive& ar, unsigned int version) const
	{
		size_t count = list->size();
		ar & count;
		for(T& t : *list)
			ar & t;
	}

	template<class Archive>
	void load(Archive& ar, unsigned int version)
	{
		size_t count;
		ar & count;
		for(size_t i=0;i<count;++i)
		{
			store->push_back(T());		  // default construct
			ar & store->back();   		  // deserialize
			list->push_back(store->back()); // store in intrusive::list
		}
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
	Store<T,Alloc>* store;
	List<T>*		list;
};


#endif /* LISTSERIALIZER_HPP_ */
