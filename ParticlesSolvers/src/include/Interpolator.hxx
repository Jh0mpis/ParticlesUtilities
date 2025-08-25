#ifndef INTERPOLATOR_HXX
#define INTERPOLATOR_HXX

#include "AMReX_Box.H"
#include "AMReX_Config.H"
#include "AMReX_MFIter.H"
#include "AMReX_REAL.H"
#include "AMReX_Random.H"
#include "AMReX_RandomEngine.H"
#include "AMReX_Scan.H"
#include "Discretizer.hxx"
#include "cctk_Types.h"
#include <array>
#include <iostream>

namespace Interpolator {

using namespace Discretize;

// #############################################################################
//                   Polynomial Interpolators
// #############################################################################

template <typename T, int Order>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
polynomial_interpolator(amrex::Array4<T const> const &f,
                        amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &x,
                        amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &plo,
                        amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &dx,
                        T &fx) {
  if constexpr (Order == 1) {
    // Linear interpolator

    // Get the coordinates in terms of the grid
    CCTK_REAL xp = (x[0] - plo[0]) / dx[0];
    CCTK_REAL yp = (x[1] - plo[1]) / dx[1];
    CCTK_REAL zp = (x[2] - plo[2]) / dx[2];

    // Compute the closest points to the position x
    // cell indexes
    int i = amrex::Math::floor(xp);
    int j = amrex::Math::floor(yp);
    int k = amrex::Math::floor(zp);

    // Compute the distances to the points
    CCTK_REAL dist_x = xp - i;
    CCTK_REAL dist_y = yp - j;
    CCTK_REAL dist_z = zp - k;

    // Compute the result
    fx = f(i, j, k) + dist_x * first_derivative<0, T, 2>(f, i, j, k, dx) +
         dist_y * first_derivative<1, T, 2>(f, i, j, k, dx) +
         dist_z * first_derivative<2, T, 2>(f, i, j, k, dx);

  } else {

    std::cout << "INFO (PARTICLESSOLVERS): Polynomial interpolation of order "
              << Order << " not yet implemented." << std::endl;
    throw std::invalid_argument(
        "Wrong order of polynomial interpolation for the first derivative");
  }
}


// #############################################################################
//                   Trilinear interpolator
// #############################################################################
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
trilinear_interpolator(amrex::Array4<T const> const &f,
                       amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &x,
                       amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &plo,
                       amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &dx,
                       T &fx) {
  // Trilinear interpolator

  // Get the coordinates in terms of the grid
  CCTK_REAL xp = (x[0] - plo[0]) / dx[0];
  CCTK_REAL yp = (x[1] - plo[1]) / dx[1];
  CCTK_REAL zp = (x[2] - plo[2]) / dx[2];

  // Compute the closest points to the position x
  // cell indexes
  int i = amrex::Math::floor(xp);
  int j = amrex::Math::floor(yp);
  int k = amrex::Math::floor(zp);

  // Compute the distances to the points
  CCTK_REAL dist_x = xp - i;
  CCTK_REAL dist_y = yp - j;
  CCTK_REAL dist_z = zp - k;
  CCTK_REAL sx[] = {CCTK_REAL(1) - dist_x, dist_x};
  CCTK_REAL sy[] = {CCTK_REAL(1) - dist_y, dist_y};
  CCTK_REAL sz[] = {CCTK_REAL(1) - dist_z, dist_z};

  // Compute the result
  fx = sx[0] * sy[0] * sz[0] * f(i, j, k) +
       sx[0] * sy[0] * sz[1] * f(i, j, k + 1) +
       sx[0] * sy[1] * sz[0] * f(i, j + 1, k) +
       sx[0] * sy[1] * sz[1] * f(i, j + 1, k + 1) +
       sx[1] * sy[0] * sz[0] * f(i + 1, j, k) +
       sx[1] * sy[0] * sz[1] * f(i + 1, j, k + 1) +
       sx[1] * sy[1] * sz[0] * f(i + 1, j + 1, k) +
       sx[1] * sy[1] * sz[1] * f(i + 1, j + 1, k + 1);
}

} // namespace Interpolator

#endif // !INTERPOLATOR_HXX
