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
    : public BaseParticleContainer<PhotonsContainer<StructType>, StructType> {

public:
  // Using BaseParticlesContainer constructor
  using Base = BaseParticleContainer<PhotonsContainer<StructType>, StructType>;
  using Base::Base;

  void evolve() override;
  void get_velocity_rhs();
  void get_position_rhs();
  void get_energy_rhs();

  void half_evolve_position();
  void evolve_velocity();
  void evolve_energy();

  void redistribute_particles();
}; // PhotonsContainer class

// ##############################################################################
//                   PhotonsContainer::METHODS DECLARATION
// ##############################################################################

template <typename StructType>
void PhotonsContainer<StructType>::get_position_rhs() {
  std::cout << this->PARTICLE_UTILITIES_INFO << "Constructing position rhs"
            << std::endl;
} // PhotonsContainer::get_position_rhs

template <typename StructType>
void PhotonsContainer<StructType>::get_velocity_rhs() {
  std::cout << this->PARTICLE_UTILITIES_INFO << "Constructing velocity rhs"
            << std::endl;
} // PhotonsContainer::get_velocity_rhs

template <typename StructType>
void PhotonsContainer<StructType>::get_energy_rhs() {
  std::cout << this->PARTICLE_UTILITIES_INFO << "Constructing energy rhs"
            << std::endl;
} // PhotonsContainer::get_energy_rhs

template <typename StructType>
void PhotonsContainer<StructType>::half_evolve_position(){
  std::cout << this->PARTICLE_UTILITIES_INFO << "x_{n+1/2} = x_n + 1/2 * dt * pos_rhs(n)"
            << std::endl;
} // PhotonsContainer::half_evolve_position

template <typename StructType>
void PhotonsContainer<StructType>::evolve_velocity(){
  std::cout << this->PARTICLE_UTILITIES_INFO <<  "v_{n+1} = RK<N>(v_n, vel_rhs)"
            << std::endl;
} // PhotonsContainer::evolve_velocity

template <typename StructType>
void PhotonsContainer<StructType>::evolve_energy(){
  std::cout << this->PARTICLE_UTILITIES_INFO << "E_{n+1} = E_n + dt * energy_rhs"
            << std::endl;
} // PhotonsContainer::evolve_energy

template <typename StructType>
void PhotonsContainer<StructType>::redistribute_particles(){
  std::cout << this->PARTICLE_UTILITIES_INFO << "Redistribute the particles to the correspondent procesor"
            << std::endl;
} // PhotonsContainer::redistribute_particles

template <typename StructType> void PhotonsContainer<StructType>::evolve() {
  std::cout << this->PARTICLE_UTILITIES_INFO << "Evolving " << this->name
            << std::endl;

  this->get_position_rhs();
  this->get_velocity_rhs();
  this->get_energy_rhs();

  this->half_evolve_position();
  this->evolve_velocity();

  this->get_position_rhs();
  this->half_evolve_position();
  this->evolve_energy();

  this->redistribute_particles();
} // PhotonsContainer::evolve

} // namespace Containers

#endif // !PARTICLESCONTAINER_H
