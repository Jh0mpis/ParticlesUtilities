# Particles Utilities

This thorn contain several utilities for the evolution of particles using the AMReX Particles container approach.

## Index

- [Particles Utilities interface](#particles-utilities-interface)
    - [Particles Container](#particles-container)
- [Files structure](#files-structure)
- [Dependencies](#dependencies)
- [Additional thorn features](#additional-thorn-features)
    - [Banned zones check](#banned-zones-check)
    - [Metric Initialization](#metric-initialization)
- [Parameters file](#parameters-file)
    - [Output parameters](#output-parameters)
    - [Banned regions parameters](#banned-regions-parameters)
    - [Metric initialization parameters](#metric-initialization-parameters)
- [Include it into your project](#include-it-into-your-project)

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
├── ParticlesContainer/
│   ├── configuration.ccl
│   ├── interface.ccl
│   ├── param.ccl
│   ├── schedule.ccl
│   └── src/
│       ├── BaseParticleContainer.hxx
│       ├── Initializers.hxx
│       ├── make.code.defn
│       ├── metrics.cxx
│       └── Particles.hxx
├── docs/
└── README.md
```

Where the  `docs` folder contains the html documentation of the project.

All the `BaseParticleContainer` and the iterator definitions are contained at the `BaseParticleContainer.hxx` file. The metrics implementations are defined inside of the `metrics.cxx` file.

## Dependencies

In order to use this thorn you need to have compiled the following thorns:

1. AMReX
1. NSIMD
1. CarpetX/Loop

and its requirements.

## Additional thorn features

The thorn also includes some extra features that can be common for all the derived particles.

### Banned zones check

Given some metric, we could find some regions were the geodesics evolution have no sense or could reach to mathematical and computing singularities. For this reason the thorn implements parameters that allows to check if the particles are crossing the banned regions. This parameters allows the user to check up to 10 different regions by specifying the $x, y, z$ coordinates and a radius. The program is going to delete the particles that enter in that 3D sphere.

This is managed by the parameters `banned_regions`, `region_<n>_position` and `region_<n>_radius` variables, where n goes from 1 to 10.

### Metric Initialization

This thorn also provides some useful implementations, such as the initialization of some space time metrics using the ADM formalism, the current implemented metrics are:

- Cartesian Minkowski.
- Schwarzschild in Isotropic coordinates.
- Schwarzschild in Schwarzschild coordinates transformed into cartesian coordinates.
- Kerr-Shild metric using cartesian coordinates.

This is managed through the parameters `metric` to define the metric type and the parameters `metric_params_i` and `metric_params_d` used for the metric parameters needed.

## Parameters file

The thorns contains the following parameters:

### Output parameters

The thorn has the following parameters for the data output:

| Parameter           | Type    | Description                                                                                          | Range                          | Default |
| ------------------- | ------- | ---------------------------------------------------------------------------------------------------- | ------------------------------ | ------- |
| particle_plot_every | INTEGER | The number of step where the simulation should print the particles data using amrex print utilities. | [0, $\infty$) 0 means no print | 0       |
| particle_tsv_every  | INTEGER | The number of step where the simulation should print the particles data using ascii.                 | [0, $\infty$) 0 means no print | 0       |

### Banned regions parameters

In order to manage the banned regions check, the thorn has defined the following parameters:

| Parameter           | Type    | Description                                        | Range                                    | Default             |
| ------------------- | ------- | -------------------------------------------------- | ---------------------------------------- | ------------------- |
| banned_regions      | INTEGER | Number of banned regions for the simulation.                         | [0, 10]                               | 0                    |
| region_<n>_position | REAL[3] | Size 3 array containing the coordinates of the banned region center. | any                                   | 0.0                  |
| region_<n>_radius   | REAL | Banned region's radius.                                                 | any                                   | 0.0                  |

### Metric initialization parameters

For the metric parameters we have the following ones:

| Parameter           | Type        | Description                                                       | Range                                                                         | Default             |
| ------------------- | ----------- | -------------------------------------------------- | ---------------------------------------- | ------------------- |
| metric              | KEYWORD     | Keyword for the type of metric the user wants to initialize.      | "Minkowski", "Isotropic Schwarzschild", "Schwarzschild", "Kerr-Shild", "none" | none |
| metric_params_d     | REAL[10]    | double type parameters needed to pass to the metric Initializer.  | any                                                                           | 0.0  |
| metric_params_i     | INTEGER[10] | integer type parameters needed to pass to the metric Initializer. | any                                                                           | 0    |

## Include it into your project

You could add this thorn to your project by adding the following lines to your EinsteinToolkit thornlist:

```bash
# PUInX
!TARGET   = $ARR
!TYPE     = git
!URL      = https://github.com/Jh0mpis/ParticlesUtilities.git
!REPO_BRANCH = dev
!REPO_PATH = $2
!CHECKOUT =
ParticlesUtilities/ParticlesContainer
```

and then you can execute the EinsteinToolkit re-build command that you usually use, for instance:

```bash
./simfactory/bin/sim build -j4 ...
```

This thorn is an abstract implementation for other thorns, so even after you code compiles, this thorn have to do nothing. Other thorns are going to be the ones that add particles to your simulation.
