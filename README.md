# dendritic-spine-geometry
A geometric computational model exploring how dendritic spine organization influences synapse formation, with parameter sweeps comparing regular and irregular dendritic architectures.

## Motivation
Abnormal dendritic spine morphology is a hallmark of several neurodevelopmental disorders and is associated with impaired synaptic connectivity. This model examines how purely geometric factors, such as spine length, orientation, spacing, and axon distance, constrain opportunities for synapse formation.

## Model Overview
- Dendrite represented as a linear segment with spines
- Axon modeled as a parallel line at varying distances
- Synapse formation requires:
  1. Spine head within a distance threshold of the axon
  2. Sufficient angular alignment toward the axon

Both regular and irregular dendritic architectures are simulated.

## Key Experiments
- Parameter sweep over axon–dendrite distance
- Parameter sweep over spine angle variability
- Comparison of regular vs. irregular spine organization

## Files
- `spine_model_stage1_1.m`: Core geometric model and visualization
- `spine_sweep.m`: Parameter sweep experiments and statistical analysis

## Future Directions
The model is designed to support future integration with biophysical neuron models
(e.g., Pinsky–Rinzel) to link dendritic geometry with electrical excitability.
