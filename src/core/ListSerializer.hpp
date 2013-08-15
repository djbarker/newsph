#ifndef LISTSERIALIZER_HPP_
#define LISTSERIALIZER_HPP_

#include <boost/serialization/serialization.hpp>
#include <memory>

/*
 * A class which allows us to serialize boost::intrusive::list objects
 * by wrapping them in an object which uses a std container as a store.
 */
template<class List, class Store, class T>
class ListSerializer
{
public:
	ListSerializer() = default;

	ListSerializer(Store& s, List& l)
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
	Store*	store;
	List*	list;
};


#endif /* LISTSERIALIZER_HPP_ */
