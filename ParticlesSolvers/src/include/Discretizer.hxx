#ifndef DISCRETIZER_HXX
#define DISCRETIZER_HXX

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <loop.hxx>
#include <stdexcept>

namespace Discretize {

template <int Dim, typename T, int Order>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
first_derivative(amrex::Array4<T const> const &gf, int i, int j, int k,
                 int comp, amrex::GpuArray<T, 3> const &dxi) {
  const unsigned int is_x = Dim == 0 ? 1 : 0;
  const unsigned int is_y = Dim == 1 ? 1 : 0;
  const unsigned int is_z = Dim == 2 ? 1 : 0;

  const T m1 = gf(i - is_x, j - is_y, k - is_z, comp);
  const T p1 = gf(i + is_x, j + is_y, k + is_z, comp);

  if constexpr (Order == 2) {

    return (-m1 + p1) / (2. * dxi[Dim]);

  } else if constexpr (Order == 4) {

    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);

    return (-8. * m1 + m2 + 8. * p1 - p2) / (12. * dxi[Dim]);

  } else if constexpr (Order == 6) {

    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);
    const T m3 = gf(i - 3 * is_x, j - 3 * is_y, k - 3 * is_z, comp);
    const T p3 = gf(i + 3 * is_x, j + 3 * is_y, k + 3 * is_z, comp);

    return (-45. * m1 + 9. * m2 - m3 + 45. * p1 - 9. * p2 + p3) /
           (60. * dxi[Dim]);

  } else if constexpr (Order == 8) {

    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);
    const T m3 = gf(i - 3 * is_x, j - 3 * is_y, k - 3 * is_z, comp);
    const T p3 = gf(i + 3 * is_x, j + 3 * is_y, k + 3 * is_z, comp);
    const T m4 = gf(i - 4 * is_x, j - 4 * is_y, k - 4 * is_z, comp);
    const T p4 = gf(i + 4 * is_x, j + 4 * is_y, k + 4 * is_z, comp);

    return (-672. * m1 + 168. * m2 - 32. * m3 + 3. * m4 + 672. * p1 -
            168. * p2 + 32. * p3 - 3. * p4) /
           (840. * dxi[Dim]);

  } else {
    std::cout << "INFO (PARTICLESSOLVERS): Discretization of the first "
                 "derivative at order "
              << Order << " not yet implemented." << std::endl;
    throw std::invalid_argument(
        "Wrong order of discretization for the first derivative");
  }
}

// #############################################################################
//                   Second derivative discretization
// #############################################################################

template <int Dim, typename T, int Order>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
second_derivative(amrex::Array4<T const> const &gf, int i, int j, int k,
                  int comp, amrex::GpuArray<T, 3> const &dxi) {

  const unsigned int is_x = Dim == 1 ? 1 : 0;
  const unsigned int is_y = Dim == 1 ? 1 : 0;
  const unsigned int is_z = Dim == 3 ? 1 : 0;

  const T m1 = gf(i - is_x, j - is_y, k - is_z);
  const T gf0 = gf(i, j, k);
  const T p1 = gf(i + is_x, j + is_y, k + is_z);

  if constexpr (Order == 2) {

    return (m1 - 2. * gf0 + p1) / (dxi[Dim] * dxi[Dim]);

  } else if constexpr (Order == 4) {

    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);

    return (-m2 + 16. * m1 - 30. * gf0 + 16. * p1 - p2) /
           (12. * dxi[Dim] * dxi[Dim]);

  } else if constexpr (Order == 6) {

    const T m3 = gf(i - 3 * is_x, j - 3 * is_y, k - 3 * is_z, comp);
    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);
    const T p3 = gf(i + 3 * is_x, j + 3 * is_y, k + 3 * is_z, comp);

    return (2. * m3 - 27. * m2 + 270. * m1 - 490. * gf0 + 270. * p1 - 27. * p2 +
            2. * p3) /
           (180. * dxi[Dim] * dxi[Dim]);

  } else if constexpr (Order == 8) {

    const T m4 = gf(i - 4 * is_x, j - 4 * is_y, k - 4 * is_z, comp);
    const T m3 = gf(i - 3 * is_x, j - 3 * is_y, k - 3 * is_z, comp);
    const T m2 = gf(i - 2 * is_x, j - 2 * is_y, k - 2 * is_z, comp);
    const T p2 = gf(i + 2 * is_x, j + 2 * is_y, k + 2 * is_z, comp);
    const T p3 = gf(i + 3 * is_x, j + 3 * is_y, k + 3 * is_z, comp);
    const T p4 = gf(i + 4 * is_x, j + 4 * is_y, k + 4 * is_z, comp);

    return (-9. * m4 + 128. * m3 - 1008. * m2 + 8064. * m1 - 14350. * gf0 +
            8064. * p1 - 1008. * p2 + 128. * m3 - 9. * m4) /
           (5040. * dxi[Dim] * dxi[Dim]);

  } else {

    std::cout << "INFO (DISCRETIZER): Discretization of the second derivative "
                 "at order "
              << Order << " not yet implemented." << std::endl;
    throw std::invalid_argument(
        "Wrong order of discretization for the second derivative");

  } // else Order not Implemented
}

} // namespace Discretize

#endif // !DISCRETIZER_HXX
