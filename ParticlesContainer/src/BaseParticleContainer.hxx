/**
 * \file BaseParticleContainer.hxx
 *
 * \brief Contains the BaseParticlesContainer and ParticleIterator classes
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
#include <cctk.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_core.h>
#include <iostream>
#include <string>

// Starting the namespace
namespace BaseContainer {

/**
 * \brief BaseParticleContainer abstract class definition.
 *
 * The BaseParticleContainer class is an abstract class that defines the methods
 * that have to be defined for the other <Particle>Containers derived classes.
 * This is templated on the new Container and the struct that defines the
 * parameters of the particles.
 */
template <typename OtherContainer, typename StructType>
class BaseParticleContainer
    : public amrex::AmrParticleContainer<0, 0, StructType::n_attributes, 0> {
public:
  // Name of the particles
  const std::string name = StructType::name;
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
   * \brief Initialize the particles given a function with a tensor product of
   * the number of particle per each cell direction..
   *
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
    initializer_function(static_cast<OtherContainer &>(*this), nppc);
  };

  /**
   * \brief Initialize the particles given a function with a fixed number of
   * particles on each cell.
   *
   * The initialize method receives a function and use it to assign the initial
   * conditions over the particles. This function can be implemented by other
   * users.
   * @param initializer_function Function that receives a BaseParticleContainer
   * instance and a vector of size AMREX_SPACEDIM with the number of particles
   * on each cell.
   * @param num_particles_per_cell Number of particles on each cell.
   * @param metric 3+1 ADM metric.
   * @param level AMR level.
   */
  template <typename Function>
  void initialize(Function initializer_function,
                  const int num_particles_per_cell,
                  const amrex::MultiFab &metric, const int &level) {
    initializer_function(static_cast<OtherContainer &>(*this),
                         num_particles_per_cell, metric, level);
  };

  /**
   * \brief Initialize the particles given a function with a fixed number of
   * particle on the container.
   *
   * The initialize method receives a function and use it to assign the initial
   * conditions over the particles. This function can be implemented by other
   * users.
   * @param initializer_function Function that receives a BaseParticleContainer
   * instance and a vector of size AMREX_SPACEDIM with the number of particles
   * on each cell.
   * @param num_particles_per_container Number of particles on the container.
   * @param metric 3+1 ADM metric.
   */
  template <typename Function>
  void initialize(Function initializer_function,
                  const int num_particles_per_container,
                  const amrex::MultiFab &metric) {
    initializer_function(static_cast<OtherContainer &>(*this),
                         num_particles_per_container, metric);
  };

  /**
   * The evolve method evolve the system given the differential equations and
   * the computed rhs.
   *
   * @param lapse ADM lapse function
   * @param shift ADM shift vector
   * @param metric ADM induced metric
   * @param curv ADM extrinsic curvature
   * @param dt Time step
   * @param lev AMR level
   */
  virtual void evolve(const amrex::MultiFab &lapse,
                      const amrex::MultiFab &shift,
                      const amrex::MultiFab &metric,
                      const amrex::MultiFab &curv, const CCTK_REAL &dt,
                      const int &lev) = 0;

  void outputParticlesAscii(CCTK_ARGUMENTS, const int &plot_every,
                            const std::string &out_dir) {

    const int it = cctkGH->cctk_iteration;
    if (plot_every > 0 && it % plot_every == 0) {
      const std::string &file_name =
          out_dir + "/" + amrex::Concatenate(this->name, it);
      CCTK_VINFO(" Writing ascii file %s", file_name.c_str());

      this->WriteAsciiFile(file_name);
    }
  };

  void outputParticlesPlot(CCTK_ARGUMENTS, const int &plot_every,
                           const std::string &out_dir) {

    const int it = cctkGH->cctk_iteration;
    if (plot_every > 0 && it % plot_every == 0) {
      const std::string file_name =
          out_dir + "/" + amrex::Concatenate("particles", it);
      CCTK_VINFO(" Writing plot file %s", file_name.c_str());

      this->WritePlotFile(file_name, "particles");
    }
  };
}; // class BaseParticlesContainer

} // namespace BaseContainer

namespace Iterator {

template <typename StructType>
class ParticleIterator
    : public amrex::ParIter<0, 0, StructType::n_attributes, 0> {
public:
  using amrex::ParIter<0, 0, StructType::n_attributes, 0>::ParIter;
  using RealVector = typename amrex::ParIter<
      0, 0, StructType::n_attributes>::ContainerType::RealVector;

  const std::array<RealVector, StructType::n_attributes> &GetAttribs() const {
    return this->GetStructOfArrays().GetRealData();
  }

  std::array<RealVector, StructType::n_attributes> &GetAttributes() {
    return this->GetStructOfArrays().GetRealData();
  }

  const RealVector &GetAttribs(int comp) const {
    return this->GetStructOfArrays().GetRealData(comp);
  }

  RealVector &GetAttributes(int comp) {
    return this->GetStructOfArrays().GetRealData(comp);
  }
}; // class ParicleIterator

} // namespace Iterator

#endif // !BASEPARTICLESCONTAINER_HXX
