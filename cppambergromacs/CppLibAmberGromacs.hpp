#ifndef CPPLIBAMBERGROMACS_HPP
#define CPPLIBAMBERGROMACS_HPP

#ifdef USE_DOUBLE_PRECISION
    using Real= double;
#else
    using Real= float;
#endif

#include "ReaderFactory.hpp"
#include "Configurations/ConfigurationBulk.hpp"
#include "Configurations/ConfigurationLipid.hpp"

#endif // CPPLIBAMBERGROMACS_HPP
