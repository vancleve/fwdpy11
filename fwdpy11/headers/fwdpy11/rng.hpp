#ifndef FWDPY11_RNG_HPP__
#define FWDPY11_RNG_HPP__

#include <fwdpp/sugar/GSLrng_t.hpp>

namespace fwdpy11
{
    /*!
      Random number generator.

      This is a std::unique_ptr wrapper to a gsl_rng * initialized
      as a Mersenne twister type (gsl_rng_mt19937).
    */
    using GSLrng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
}

#endif
