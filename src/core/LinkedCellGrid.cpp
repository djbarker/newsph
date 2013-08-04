#include "LinkedCellGrid.hpp"

using namespace std;

namespace sim
{

/*
 * Really we would want to partially specialize on Dim (and possibly Padding)
 * but we cannot partially specialize members without partially specializing
 * the whole class so we are forced to use ifs. This removes some compile-
 * time safety but it's more acceptable than the template work-arounds. It
 * would also be safer to have loc as a template parameter so we could tell
 * at compile time if we were accidently calling say getPadding<Front> with
 * Dims==2 which would be invalid.
 */

template<size_t Dim, class T, size_t Padding>
list<T> LinkedCellGrid<Dim,T,Padding>::getPadding(size_t loc)
{
	if(Dim==2)
	{
		// TODO: padding
		switch(loc)
		{
		case Top:
			break;
		case Bottom:
			break;
		case Left:
			break;
		case Right:
			break;
		case Top|Left:
			break;
		case Top|Right:
			break;
		case Bottom|Left:
			break;
		case Bottom|Right:
			break;
		default:
			throw logic_error("invald padding location for 2D LinkedCellGrid!"); // see above
		}
	}
	else if(Dim==3)
	{

	}

}

} /* namespace sim */
