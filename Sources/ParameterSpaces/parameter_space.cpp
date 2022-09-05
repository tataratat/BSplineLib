#ifndef SOURCES_PARAMETERSPACE_PARAMETER_SPACE_CPP_
#define SOURCES_PARAMETERSPACE_PARAMETER_SPACE_CPP_

#include "Sources/ParameterSpaces/parameter_space.hpp"

#ifdef SPLINELIB_EXPLICIT_STATIC_LIBRARY
template class ParameterSpace<1>;
template class ParameterSpace<2>;
template class ParameterSpace<3>;
template class ParameterSpace<4>;
template class ParameterSpace<5>;
template class ParameterSpace<6>;
template class ParameterSpace<7>;
template class ParameterSpace<8>;
template class ParameterSpace<9>;
template class ParameterSpace<10>;
#endif

#endif
