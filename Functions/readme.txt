There are three main scripts in this folder:

- The main one is "full forward equations v3.r". It contains two functions that are used to run the moment-closure forward equations. The purpose of the first function is to make the main forward-equation code be as efficient memory-wise and flops-wise as possible by calculating in advance the correct indices for all the updates and to not store duplicates of, e.g., sxx(s_i, s_j) = sxx(s_j, s_i).

- The "mcmc functions -- emulator.r" is for the emulator-setup in the simulation studies.

- "custom gillespie v4.r" is used to simulate the data and is simply a custom Gillespie stochastic simulator. The version with "-logdens.r" is meant to easily accommodate a spatially varying beta.