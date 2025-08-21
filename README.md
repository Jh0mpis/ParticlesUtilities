# Particles Utilities

This thorn contain several utilities for the evolution of particles using the AMReX Particles container approach.

## Index

-

## Particles Utilities interface

The `ParticlesUtilities` thorn is an interface thorn that contains the common classes and functions for differential equation associated to systems of particles. This thorn has the following utilities implemented:

### Particles Container

The particles container thorn contains the definition of the `BaseParticleContainer` templated class. This abstract class defines the basic methods that have to be defined on each of the derived classes, this class extends from the `amrex::AmrParticleContainer`.

The `BaseParticleContainer` class is templated by the new `<Particle>Container` derived class first, and a struct object defining the particle position, velocity or momentum and another extra quantities needed by the particles. The thorn also contains some particles defined by default inside of the `Particles.hxx` header file.

This class also implements some initializers functions, depending on the initial conditions that we want to apply over the particles, the function that is going to be called to initialize the system can be also defined by the user following the structure defined inside of the `Initializers.hxx` header file.

> [!NOTE]
> The header file `TestParticlesContainer.hxx` contains an example of a derived BaseParticleContainer class. This could be removed on the release version of the thorn and may be included in the documentation. The same is the case of the `test.cxx`, this is just a file to check how could be the flow of a common application.

### Particles Solvers

The `ParticlesSolvers` contains the discretization of the first and second partial derivatives of a given function. The discretizator functions are templated depending on the direction of the derivative, the type of the function (scalar, vectorial or tensorial) and the order of the discretization (2, 4, 6 or 8).

> [!NOTE]
> This thorn should include also the interpolators and solvers depending on the needings we define.

## Files structure

The project's files are distributed in the following way:

```bash
ParticlesUtilities/
├── .gitignore
├── ParticlesContainer/
│   ├── configuration.ccl
│   ├── doc/
│   ├── interface.ccl
│   ├── par/
│   │   └── test.par
│   ├── param.ccl
│   ├── schedule.ccl
│   └── src/
│       ├── include/
│       │   ├── BaseParticleContainer.hxx
│       │   ├── Initializers.hxx
│       │   ├── make.code.defn
│       │   ├── Particles.hxx
│       │   └── TestParticlesContainer.hxx
│       ├── make.code.defn
│       └── test.cxx
├── ParticlesSolvers/
│   ├── configuration.ccl
│   ├── doc/
│   ├── interface.ccl
│   ├── param.ccl
│   ├── schedule.ccl
│   └── src/
│       ├── include/
│       │   ├── Discretizer.hxx
│       │   └── make.code.defn
│       └── make.code.defn
└── README.md
```

<!-- ## How to add it to a project -->

