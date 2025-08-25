/**
 * @file Particles.hxx
 * 
 * This file contains some predefined structs that defines the parameters
 * needed for a particle in order to be used on a <Particles>Container instance.
 */

#ifndef PARTICLES_HXX
#define PARTICLES_HXX

namespace Particles {
// Struct to managing null particles
struct BaseParticleData {
  static constexpr const char *name = "Basic Particles";
  enum {
    // Particle position
    x = 0, y, z,
    // Particle Velocity
    vx, vy, vz,
    // Number of attributes on the struct
    n_attributes,
  }; // enum
}; // struct PhotonsData

// Struct to managing null particles
struct PhotonsData {
  static constexpr const char *name = "Photons";
  enum {
    // Photon position
    x = 0, y, z,
    // Photon Velocity
    vx, vy, vz,
    // Photon Energy
    E,
    // Number of attributes on the struct
    n_attributes,
  }; // enum
}; // struct PhotonsData

// Struct to managing massive particles
struct MassiveParticleData {
  static constexpr const char *name = "Massive Particles";
  enum {
    // Massive particle position
    x = 0, y, z,
    // Massive particle velocity
    vx, vy, vz,
    // Massive particle Energy
    E,
    // Massive particle mass
    m,
    // Number of attributes on the struct
    n_attributes,
  }; // enum
}; // struct MassiveParticleData

} // namespace Particles
#endif // !PARTICLES_HXX
