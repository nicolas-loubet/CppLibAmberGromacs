#ifndef CPPLIBAMBERGROMACS_HPP
#define CPPLIBAMBERGROMACS_HPP

#include <string>

#ifdef USE_DOUBLE_PRECISION
    using Real= double;
    double RealParser(const std::string& s) { return std::stod(s); }
#else
    using Real= float;
    float RealParser(const std::string& s) { return std::stof(s); }
#endif

#include "ReaderFactory.hpp"
#include "Configurations/ConfigurationBulk.hpp"
#include "Configurations/ConfigurationLipid.hpp"

#endif // CPPLIBAMBERGROMACS_HPP
