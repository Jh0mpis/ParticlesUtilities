/**
 * @file BaseParticlesContainer.hxx
 *
 * @brief Contains the BaseParticlesContainer and ParticleIterator classes
 * along with the BaseContainer and Iterator namespaces definition.
 *
 * This file contains the BaseParticlesContainer abstract class. This
 * class extends from amrex::AmrParticleContainer and contains the particles
 * data using an array of structs (AoS).
 *
 * Also, includes the definition of the ParticleIterator class, that could be
 * helpful for the manipulation of the particles.
 */
#ifndef BASEPARTICLESCONTAINER_HXX
#define BASEPARTICLESCONTAINER_HXX

// Include libraries
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <cctk_Arguments.h>

// Starting the namespace
namespace BaseContainer {

/**
 * The BaseParticleContainer class is an abstract class that defines the methods
 * that have to be defined for the other <Particle>Containers derived classes.
 * This is templated on the new Container and the struct that defines the
 * parameters of the particles.
 */
template <typename OtherContainer, typename StructType>
class BaseParticleContainer
    : public amrex::AmrParticleContainer<0, 0, StructType::n_attributes, 0> {
public:
  // Derived class type for future casts
  using ThisContainer = OtherContainer;

  // INFO message for printing
  static constexpr const char *PARTICLE_UTILITIES_INFO =
      "INFO (ParticlesUtilities): ";
  // Name of the particles
  static constexpr const char *name = StructType::name;
  // Number of attributes per each particle
  static constexpr int n_attributes = StructType::n_attributes;

  /**
   * The BaseParticleContainer constructor initialize the AmrCore superclass
   * receiving an external AmrCore instance.
   * @param amr_core A pointer to an amrex::AmrCore instance.
   */
  BaseParticleContainer(amrex::AmrCore *amr_core)
      : amrex::AmrParticleContainer<0, 0, StructType::n_attributes, 0>(
            amr_core) {}

  /**
   * The initialize method receives a function and use it to assign the initial
   * conditions over the particles. This function can be implemented by other
   * users.
   * @param initializer_function Function that receives a BaseParticleContainer
   * instance and a vector of size AMREX_SPACEDIM with the number of particles
   * per direction on each cell.
   * @param nppc Number of particles per direction on each cell.
   */
  template <typename Function>
  void initialize(Function initializer_function,
                  const std::array<int, AMREX_SPACEDIM> nppc) {
    initializer_function(static_cast<ThisContainer &>(*this), nppc);
  };

  /**
   * The evolve method evolve the system given the differential equations and
   * the computed rhs.
   */
  virtual void evolve() = 0;
  // virtual void computeRHS() = 0;

  void outputParticlesAscii(CCTK_ARGUMENTS);
  void outputParticlesPlot(CCTK_ARGUMENTS);
}; // class BaseParticlesContainer

} // namespace BaseContainer

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
