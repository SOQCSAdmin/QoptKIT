# Quantum Optical KIT for QISKIT (QoptKIT) 
** Work in progress. Under testing. **

Quantum Optical KIT (QoptKit) library adds to QISKIT the capability of simulating and compiling ideal optical circuis. Optical circuits are defined from its ideal components connected between them. The library also provides support to establish post-selection conditions on the circuit. The aim of this library is to provide a first introduction to quantum optical technologies to QISKIT users.

**QoptKIT capabilities**:
- Create an ideal circuit made of phase shifters and beamsplitters
- Create a simulator to run it
- Run the simulation
- Translate a QoKITt outcome into qubit encoding.
- Print the output statistics. 

 
# 1. Requirements #

* Linux or MaxOsX operating system
* C++ Compiler
* GNU Make
* [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* Python3
* qiskit
* matplotlib
* numpy


# 2. Installation #
pip install git+https://github.com/SOQCSADmin/QoptKIT

This library has various C++ files that have to be compiled. It is necedary to have installed a C++ compiler and the tools Make and ar.  These are standard tools avaiable
for any Linux distribution. Use you preferred package installer to have them available in your system. QoptKit also require the C++ library [Eigen 3] to be installed in your system.
This library is also available in most of the Linux distributions.

Python libraries qiskit, matplotlib and numpy are required but they will be installed acutomatically. We recommend the use of conda but its not required.


# 3 Documentation
The documentation can be found in https://soqcsadmin.github.io/QoptKIT/

# 4. Authorship #
<b>Javier Osca</b> <br>
javier.oscacotarelo@mu.ie

<b>Jiri Vala</b> <br>
jiri.vala@mu.ie

# 5. License and Copyright #
Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. 

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the [License](./LICENSE.TXT). 

