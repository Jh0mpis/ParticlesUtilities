# Particles Utilities

This thorn contain several utilities for the evolution of particles using the AMReX Particles container approach.

## Index

- [Particles Utilities interface](#particles-utilities-interface)
    - [Particles Container](#particles-container)
- [Files structure](#files-structure)

## Particles Utilities interface

The `ParticlesUtilities` thorn is an interface thorn that contains the common classes and functions for differential equation associated to systems of particles. This thorn has the following utilities implemented:

### Particles Container

The particles container thorn contains the definition of the `BaseParticleContainer` templated class. This abstract class defines the basic methods that have to be defined on each of the derived classes, this class extends from the `amrex::AmrParticleContainer`.

The `BaseParticleContainer` class is templated by the new `<Particle>Container` derived class first, and a struct object defining the particle position, velocity or momentum and another additional quantities needed by the particles. The thorn also contains some particles defined by default inside of the `Particles.hxx` header file this file can be used as an example file for new particles definitions.

This class also implements some initializers functions, depending on the initial conditions that we want to apply over the particles, the function that is going to be called to initialize the system can be also defined by the user following the structure defined inside of the `Initializers.hxx` header file, again this file can be used as a starting point for custom initializer functions.

> [!NOTE]
> The functions defined inside `Initializers.hxx` assumes Particle structs similar to the ones on `Particles.hxx`, some of it could not work if you have different variables or names on your particle struct.

## Files structure

The project's files are distributed in the following way:

```bash
ParticlesUtilities/
├── ParticlesContainer
│   ├── configuration.ccl
│   ├── interface.ccl
│   ├── param.ccl
│   ├── schedule.ccl
│   └── src
│       ├── BaseParticleContainer.hxx
│       ├── Initializers.hxx
│       ├── make.code.defn
│       └── Particles.hxx
├── docs/
└── README.md
```

<!-- ## How to add it to a project -->

