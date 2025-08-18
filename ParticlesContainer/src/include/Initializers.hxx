#include "AMReX_Config.H"
#include <iostream>

template <typename ParticleContainerClass>
void test_initializer(ParticleContainerClass &pc, const int nppc){
  std::cout << pc.PARTICLE_UTILITIES_INFO << "Initializing "<<pc.name<< " with "
            << pc.n_attributes << " attributes." << std::endl;
  std::cout << pc.PARTICLE_UTILITIES_INFO << "Initializing Photons in "
            << AMREX_SPACEDIM << " dimensions." << std::endl;
}
