#if !defined(TIMESTEP_H)
#define TIMESTEP_H

#include <cstddef>

namespace timestep {

  template<typename Real>
    void timestep(size_t n_asteroids,
                  Real mu, Real dt,
                  Real *__restrict__ x,
                  Real *__restrict__ y,
                  Real *__restrict__ z,
                  Real *__restrict__ u,
                  Real *__restrict__ v,
                  Real *__restrict__ w);

  template<typename Real>
    void position_timestep(size_t n_asteroids,
                           Real dt,
                           Real *__restrict__ x,
                           Real *__restrict__ y,
                           Real *__restrict__ z,
                           Real *__restrict__ u,
                           Real *__restrict__ v,
                           Real *__restrict__ w);
}

#endif
