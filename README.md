This repository contains the code and data to reproduce the results of the Mascot-Skyline manuscript.
The code for MASCOT-Skyline is part of the MASCOT package available through the BEAST2 package manager. The source code for MASCOT
can be found at github.com/nicfel/MASCOT and a tutorial can be found at github.com/nicfel/MASCOT-tutorial.

# MASCOT-Skyline integrates population and migration dynamics to enhance phylogeographic reconstructions

Nicola F. Müller^a,b,1, Remco R. Bouckaert^c, Chieh-Hsi Wu^d, Trevor Bedford^b,e

- ^a Division of HIV, ID and Global Medicine, University of California San Francisco, San Francisco, USA
- ^b Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, USA
- ^c Centre for Computational Evolution, The University of Auckland, New Zealand
- ^d School of Mathematical Sciences, University of Southampton, UK
- ^e Howard Hughes Medical Institute, Seattle, USA

*^1 Corresponding author*

**Contact:** <nicola.felix.mueller@gmail.com>

**Abstract:** Phylodynamic methods can quantify temporal and spatial transmission dynamics of infectious diseases from information contained in phylogenetic trees. Usually, phylodynamic methods infer spatial or temporal transmission dynamics separately, leading to biased inferences and limiting their application to study disease spread. Here, we introduce a structured coalescent skyline approach, MASCOT-Skyline, to quantify spatial transmission patterns of infectious diseases and how population sizes and migration rates change over time. We model the effective population size dynamics in different locations using a non-parametric function, allowing us to approximate a range of population size dynamics. We implemented the inference of non-parametric population size dynamics as part of the Bayesian phylodynamics platform BEAST2 and the software package MASCOT. Using a range of data sets and simulations, we show that both temporal and spatial dynamics should be modeled to provide accurate inferences, even when only one or the other is of interest. Current methods that model either spatial or temporal transmission dynamics, but not both simultaneously, are biased in various situations. However, accounting for both simultaneously, we can retrieve complex temporal dynamics across different locations from pathogen genome data while providing accurate estimates of the transmission rates between those locations.