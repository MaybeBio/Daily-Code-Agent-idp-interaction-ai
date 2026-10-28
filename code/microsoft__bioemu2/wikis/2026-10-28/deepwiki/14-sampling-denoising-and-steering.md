# Sampling, Denoising and Steering

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/propagation.md](docs/propagation.md)
- [docs/sampling.md](docs/sampling.md)
- [docs/steering.md](docs/steering.md)

</details>



The sampling, denoising, and steering subsystem of BioEmu2 coordinates the transformation of noise and conditioning inputs into physical biomolecular structures. It encompasses Elucidating Diffusion Models (EDM) samplers, task-specific operational modes, inference-time physical and geometric steering potentials, spatial contact conditioning, and implicit transfer operator time propagation. 

This is a parent page providing a high-level overview of the denoising stack. For deep technical implementations, see the child pages:
- [Denoising Samplers and EDM](#4.1)
- [Inference-Time Steering and Potentials](#4.2)
- [Contact Conditioning](#4.3)
- [Time Propagation (Implicit Transfer Operators)](#4.4)

```mermaid
graph TD
    A["NaturalLanguage:EquilibriumSampling"] --> B["CodeEntity:EmuSession.sample"]
    B --> C["CodeEntity:EDMSampler"]
    C --> D["CodeEntity:ScoreModel"]
    
    E["NaturalLanguage:StructurePrediction"] --> F["CodeEntity:EmuSession.mode='structure'"]
    F --> C

    Sources: [docs/sampling.md:107-120]()
```

## Denoising Samplers and EDM

BioEmu2 frames generation as a stochastic diffusion process over atomic coordinates governed by score-based models [docs/sampling.md:147-151](). The foundational sampling loops are implemented around Elucidating Diffusion Models (EDM) principles. Depending on user intent, the system operates in either **equilibrium** mode or **structure prediction** mode. Equilibrium mode targets an unbiased representation of the Boltzmann ensemble, whereas structure prediction applies heuristic sharpening to produce single, highly folded reference geometries [docs/sampling.md:107-120]().

To handle rotational and translational invariance during the reverse diffusion trajectory, the denoising engine integrates specialized modules for augmentation, rigid body alignment, and context-parallel synchronization. 

For full details on the underlying samplers, noise schedules, alignment routines, and mode configurations, see [Denoising Samplers and EDM](#4.1).

Sources: [docs/sampling.md:107-151](), [docs/propagation.md:163-165]()

## Inference-Time Steering and Potentials

To influence generated conformational ensembles without retraining the underlying model, BioEmu2 provides an inference-time steering subpackage [docs/steering.md:16-22](). Steering modifies the sampling trajectory by tilting the base model distribution toward a user-defined target distribution using scalar energy potentials $U(x)$ [docs/steering.md:34-42](). 

The package supports gradient-based guidance and Feynman-Kac particle resampling via specialized sampler classes such as `GuidedFKSampler` and `FKCSampler` [docs/steering.md:78-81](). Default physicality potentials automatically penalize atomic clashes, incorrect bond lengths, and distorted geometry, ensuring chemically sensible outputs even for large macromolecular systems [docs/steering.md:61-66]().

```mermaid
graph TD
    A["NaturalLanguage:PhysicalitySteering"] --> B["CodeEntity:GuidedFKSampler"]
    A2["NaturalLanguage:ParticleResampling"] --> C["CodeEntity:FKCSampler"]
    B --> D["CodeEntity:Potential"]
    C --> D
    D --> E["CodeEntity:BondLengthPotential"]
    D --> F["CodeEntity:ClashPotential"]

    Sources: [docs/steering.md:78-96]()
```

For comprehensive documentation on custom potentials, RDKit constraints, and sampler configurations, see [Inference-Time Steering and Potentials](#4.2).

Sources: [docs/steering.md:16-96]()

## Contact Conditioning

Spatial constraints and contact maps can be passed directly into the generation pipeline to restrict or guide the layout of polymer chains [docs/sampling.md:143-148](). Unlike post-hoc steering potentials, contact conditioning injects distance or pairing constraints directly into the feature representation and trunk conditioning layers prior to diffusion sampling [docs/steering.md:28-33]().

This mechanism enables users to enforce experimentally derived interaction restraints (such as cross-linking mass spectrometry or mutagenesis data) while allowing the emulator to generate physically relaxed tertiary structures around those constraints.

For further information on contact map representations, conditioning APIs, and constraint examples, see [Contact Conditioning](#4.3).

Sources: [docs/sampling.md:143-148](), [docs/steering.md:28-33]()

## Time Propagation (Implicit Transfer Operators)

Beyond unconditioned generation and equilibrium sampling, BioEmu2 supports time propagation using implicit transfer operators. Given a starting conformation $x_0$, the model can sample future endpoints $x_\tau$ at a specified lag time $\tau$ (`lag_time_ns`), effectively simulating molecular dynamics over extended timescales without stepping through continuous short-time integration paths [docs/propagation.md:16-28]().

The propagation framework supports single-endpoint sampling per frame, multi-endpoint branching, and time-ordered trajectory construction via `simulate_trajectory()` [docs/propagation.md:16-20, 150-162](). Reusing trunk embeddings across sequential propagation steps ensures computational efficiency during extended temporal simulations [docs/propagation.md:26-28].

For detailed guides on lag-time configuration, external coordinate mapping, and trajectory generation workflows, see [Time Propagation (Implicit Transfer Operators)](#4.4).

Sources: [docs/propagation.md:16-28, 150-162]()

---