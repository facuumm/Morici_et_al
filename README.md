# Morici_et_al

for Dorsal-Ventral HPC Coordination

These codes have been developed to study dorsal-ventral  hippocampal coordination at different levels (single units, assemblies,  LFP).

In the "toolbox" folder, you'll find both general and specific functions.

In the main folder, various pipelines have been provided to perform the analysis.

The codebase is created using MATLAB, and comprehensive  documentation, along with line-by-line descriptions, has been provided  for each section of the code.

This project is a work in progress, so please feel free to contact me if you discover any bugs.

[facundo.morici@inserm.fr](mailto:facundo.morici@inserm.fr) / [faq.morici@gmail.com](mailto:faq.morici@gmail.com)







**---------------------------------------------------------------**

**List of functions**

**---------------------------------------------------------------**

**General**

- waveform_parameters*: Extract waveform parameters for Pyr/Int classification. inside ‘tools_for_SU’

- spike_train_construction*: Spike Trains matrix construction. inside ‘tools_for_SU’

- ShuffleSpks: shuffle spikes 100 keeping the inter-spike intervals. inside ‘tools_for_SU’

- merge_events*: merge events that occur close in time. inside ‘tools_for_LFP’

- inner*: inner product of two vectors. inside ‘tools_for_assemblies’

- *gaussfilt*: Gaussian filter for time series data. inside ‘tools_for_assemblies’

- count_spks*: count spikes occurring within an event. inside ‘tools_for_SU’

**Analysis**

- *ThetaModulation*: This function determines if SUs are phase-locked to theta rhythm.

- *RipplesModulation*: Is useful to determine if the SUs are modulated by ripples using a Poisson Test. inside ‘tools_for_SU’

- *SU_responsivness*: Firing Rate tuning curve calculation. Locked to an event. inside ‘tools_for_SU’

- *spatial_info_rate*: Spatial Information calculation based on Skaags et al 1993. inside ‘tools_for_SU’

- *sparsity_info*: Sparsity information calculation based on Skaags et al 1993. inside ‘tools_for_SU’

- *SimilarityIndex*: Similarity Index between assemblies patterns based on Almeida-Filho et al 2014. inside ‘tools_for_assemblies’

- *reactivation_strength*: Reactivation Strength calculation based on van de ven et al (2016). inside ‘tools_for_assemblies’

- *FiringMap_LinearTrack*: Firing curve tuned to space. inside ‘tools_for_SU’

- assembly_patternsJFM: assemblies detection based on Lopes-dos-Santos algorithm. Modified to determine assemblies members depending on different thresholding methods. inside ‘tools_for_assemblies’

  **general_for_analyses**

  - *Threshold_xcorr*: Assign random positions to time2 vector and construct a ccg using the times1 vector. inside ‘general_tools’
  - *Threshold_CCG_Random*: Similar to Threshold_xcorr.m. inside ‘general_tools’
  - *SkaggsRandom*: 90-quantile definition from a Skaags Surrogate distribution. inside ‘tools_for_SU’
  - *SkaggsRandomFMT*: Similar to *SkaggsRandom* but using FMAToolbox. inside ‘tools_for_SU’
  - *SimilaritySurrogate*: This function determines the 99-th percentile from a surrogate Similarity Index distribution using different assemblies patterns. inside ‘tools_for_assemblies’
  - *marchenko*: Marchenko Pastur distribution using SpikeTrains. inside ‘tools_for_assemblies’
  - *assembly_recruitment*: percentage of events recruited by different assemblies. inside ‘tools_for_assemblies’
  - *Within_pc*: Calculates Spatial correlation, Firing Rate Change, and Rate overlap within the same session, controlling the size of the sample for each group. inside ‘tools_for_SU’

**Plotting**

- *triggered_activity*: Reactivation Strength calculation of different assemblies surrounding an event for plotting. inside ‘tools_for_plotting’

- *Ripple_PHIST*: Tuning curve of SUs to ripples. inside ‘tools_for_plotting’

- PHIST*: Similar to *Ripple_PHIST*. inside ‘tools_for_plotting’

- *RasterPlot*: Raster plot construction using two different time series. inside ‘tools_for_plotting’

**Decoding**

- *bayesian_replay*: Probability calculation of being at different positions using spiking activity based on Bayesian probability. inside ‘tools_for_decoder’

  **general_for_decoding**

  - *FindReplay*: Find putative-replay events using MU activity. inside ‘tools_for_decoder’