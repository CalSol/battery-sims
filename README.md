# Battery simulations

Run `setup` once per MATLAB session from the repository root. It adds the source,
Simscape model assets, prototypes, and diagnostics to the MATLAB path without
changing the persistent user path.

## Layout

- `src/analyses` — supported analytical simulation entry points:
  `batterySimFramework_36s9p` and `thermalRunawayStudy_36s9p`.
- `src/builders` — Simscape model builders, including
  `packElectricalModel_36s9p`.
- `src/visualization` — circuit and pack-layout visualization builder.
- `models/simscape/assemblies` — Simscape libraries, generated package assets,
  and their associated MAT files, grouped by assembly.
- `prototypes` — exploratory thermal and PDE scripts; these are not part of the
  supported analysis interface.
- `notebooks` — MATLAB live scripts retained as exploratory records.
- `tools/simscape-probes` — one-off Simscape library and port diagnostics.
- `sim_outputs` — existing simulation outputs and exported reference figures.

`setup.m` is the only path bootstrap file. New supported MATLAB entry points
belong under `src`; generated models and package assets belong under
`models/simscape/assemblies`; throwaway investigations belong under
`prototypes` or `tools`.

## Design boundary

The repository currently uses functions and scripts rather than classes. No
object-oriented wrapper layer was introduced: it would add ceremony without
changing the simulations. If shared state later warrants OOP, keep it to small
configuration/value classes and independent services rather than deep class
hierarchies. This preserves a scalable public surface while keeping the
numerical models straightforward to inspect.
