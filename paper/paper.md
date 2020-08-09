---
title: 'Excimontec v1.0: An Open-Source Software Tool for Kinetic Monte Carlo Simulations of Organic Electronic Devices'
tags:
  - kinetic Monte Carlo
  - C++
  - organic photovoltaics
  - organic semiconductors
  - exciton diffusion
  - charge recombination
  - charge transport
  
authors:
 - name: Michael C. Heiber
   orcid: 0000-0002-1567-5663
   affiliation: "1"
affiliations:
 - name: Center for Hierarchical Materials Design (CHiMaD), Northwestern University, Evanston, Illinois 60208, USA
   index: 1
date: 4 May 2020
bibliography: paper.bib
---

# Summary

For over three decades, kinetic Monte Carlo (KMC) simulations have been a powerful computational tool to help understand and optimize organic semiconductor devices, especially photovoltaics, light-emitting diodes, transistors, and thermoelectrics.[@baranovskii2014pssb; @groves2017rpp; @heiber2019chapter; @zuo2019aem]
KMC simulations use mechanistic models for how excitons and polarons are created, migrate through, and are then eventually removed from the semiconductor layer of a device and can capture the complex interactions between performance and spatial structure that is often not possible using continuum drift-diffusion models. 
This can then be used to probe a wide variety of phenomena in organic electronics devices, including exciton diffusion and quenching, charge transport, and charge recombination at the full device scale while retaining details regarding nanoscale inhomogeneities. 
Despite the clear utility of the method, no widespread or standardized software tools have taken hold in the community. 
Instead, many research groups around the world have maintained private codebases of varying complexity, efficiency, and reliability. 
As a result, there have been large barriers to entry for new researchers and a lot of repeated effort throughout the community that would have been much better off applied to pushing the capabilities of the technique and further refining the physical models.

``Excimontec`` is designed to be a well-tested, optimized, reliable, and accessible open-source tool for performing KMC simulations of organic electronic devices. 
v1.0 has a particular focus on organic photovoltaic device modeling and can utilize complex bulk heterojunction morphologies generated using the ``Ising_OPV`` tool.[@heiber2018joss]
v1.0 comes with five different simulation tests: exciton diffusion, time-of-flight charge transport, internal quantum efficient, dynamics, and steady state charge transport. 
See the user manual for a more in-depth description each simulation test, including some examples of what each test can be used for. 
The software has been developed in modern C++ and is optimized for efficient execution on high performance computing clusters using MPI. 
This software package uses object-oriented design and extends the ``KMC_Lattice`` framework.[@heiber2019joss] 
The code includes rigorous unit and validation testing with ``googletest``, continuous integration testing with ``TravisCI``, and API documentation generated using ``Doxygen``. 
The source code for ``Excimontec v1.0`` is archived with Zenodo.[@heiber2020excimontec1.0.0]


# Acknowledgments

This work was developed under the financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).  Thank you to Dr. Dean M. DeLongchamp for providing access to NIST's Raritan computing cluster, which was helpful with software development and testing.

# References
