#include "WendlandQuintic.hpp"

namespace sim
{
namespace kernels
{

template<> const quantity<number> WendlandQuintic<2>::C(0.557042300821633675191093171803800267120608760091597570616835); // 7/(4pi)
template<> const quantity<number> WendlandQuintic<3>::C(0.417781725616225256393319878852850200340456570068698177962626); // 21/(16pi)

}
}
