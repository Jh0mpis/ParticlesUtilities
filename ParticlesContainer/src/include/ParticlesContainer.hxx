/**
 * @file ParticlesContainer.hxx
 *
 * @brief Contains the Iterator and Particles namespaces.
 *
 * This file contains the definition of the Iterator and Particles namespaces.
 *
 * The Iterator namespace contains the definition to a convenient class to
 * iterate over the particles: the ParticleIterator class. The Particles
 * namespace defines the different <Particles>Container classes, where
 * <Particle> depend on the equations to solve. This is the concrete class that
 * implements the abstract methods defined on the abstract
 * BaseParticlesContainer class.
 */

#ifndef PARTICLESCONTAINER_H
#define PARTICLESCONTAINER_H

// Import libraries
#include "BaseParticlesContainer.hxx"
#include <AMReX_AmrParticles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <cctk_Arguments.h>

namespace Containers {

// #############################################################################
//                   PhotonsContainer::CLASS INITIALIZATION
// #############################################################################
using namespace BaseContainer;

template <typename StructType>
class PhotonsContainer
    : public BaseParticleContainer<PhotonsContainer<StructType>,StructType> {

public:
    // Using BaseParticlesContainer constructor
  using Base = BaseParticleContainer<PhotonsContainer<StructType>, StructType>;
  using Base::Base;

  void evolve() override;
}; // PhotonsContainer class

// ##############################################################################
//                   PhotonsContainer::METHODS DECLARATION
// ##############################################################################

template <typename StructType> void PhotonsContainer<StructType>::evolve() {
  std::cout << this->PARTICLE_UTILITIES_INFO << "Evolving "<<this->name
            << std::endl;
} // PhotonsContainer::evolve

} // namespace Containers

#endif // !PARTICLESCONTAINER_H
