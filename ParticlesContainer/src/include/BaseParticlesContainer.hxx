/**
 * @file BaseParticlesContainer.hxx
 *
 * @ brief Contains the BaseParticlesContainer class definition.
 *
 * This file contains the BaseParticlesContainer abstract class. This
 * class extends from amrex::AmrParticleContainer and contains the particles
 * data using an array of structs (AoS).
 */
#ifndef BASEPARTICLESCONTAINER_HXX
#define BASEPARTICLESCONTAINER_HXX

// Include libraries
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <cctk_Arguments.h>

// Starting the namespace
namespace BaseContainer {

template <typename OtherContainer, typename StructType>
class BaseParticleContainer
    : public amrex::AmrParticleContainer<0, 0, StructType::n_attributes, 0> {
public:
using ThisContainer = OtherContainer;
  static constexpr const char *PARTICLE_UTILITIES_INFO =
      "INFO (ParticlesUtilities): ";
  static constexpr const char *name = StructType::name;
  static constexpr int n_attributes = StructType::n_attributes;
  BaseParticleContainer(amrex::AmrCore *amr_core)
      : amrex::AmrParticleContainer<0, 0, StructType::n_attributes, 0>(
            amr_core) {}

  template <typename Function>
  void initialize(Function initializer_function, const int nppc) {
    initializer_function(static_cast<ThisContainer&>(*this), nppc);
  };

  virtual void evolve() = 0;
  // virtual void computeRHS() = 0;

  void outputParticlesAscii(CCTK_ARGUMENTS);
  void outputParticlesPlot(CCTK_ARGUMENTS);
}; // class BaseParticlesContainer

} // namespace BaseParticlesContainer

namespace Iterator {

template <typename StructType>
class ParticleIterator
    : public amrex::ParIter<0, 0, StructType::n_attributes, 0> {
public:
  using Base = amrex::ParIter<0, 0, StructType::n_attributes, 0>;
  using Base::ParIter;
  using RealVector = typename amrex::ParIter<
      0, 0, StructType::n_attributes>::ContainerType::RealVector;

  const std::array<RealVector, StructType::n_attributes> &GetAttribs() const {
    return this->GetStructOfArrays().GetRealData();
  }

  std::array<RealVector, StructType::n_attributes> &GetAttribs() {
    return this->GetStructOfArrays().GetRealData();
  }

  const RealVector &GetAttribs(int comp) const {
    return this->GetStructOfArrays().GetRealData(comp);
  }

  RealVector &GetAttribs(int comp) {
    return this->GetStructOfArrays().GetRealData(comp);
  }
}; // class ParicleIterator

} // namespace Iterator

#endif // !BASEPARTICLESCONTAINER_HXX
