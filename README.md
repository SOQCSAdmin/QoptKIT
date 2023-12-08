# Quantum Optical KIT for QISKIT (QoptKIT) 

<p align="justify" >Quantum Optical KIT (QoptKit) library adds to QISKIT the capability of simulating and compiling quantum optical circuits. Optical circuits are defined by interconnected ideal components. The library also provides support to establish conditions for post-selection which is required in optical circuits. The aim of this library is to introduce quantum optical technologies to QISKIT users.</p>

**QoptKIT capabilities**:

- Create a circuit made of phase shifters and beamsplitters
- Create a simulator 
- Run the simulation
- Translate a QoptKIT outcome into qubit encoding.
- Print the output statistics. 

 
# 1. Requirements #

* Linux or MacOsX operating system
* C++ Compiler
* GNU Make
* [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* Python3
* qiskit


# 2. Installation #
pip install git+https://github.com/SOQCSADmin/QoptKIT

<p align="justify">This library has various C++ files that have to be compiled. It is necessary to have installed a C++ compiler and the tools make and ar.  These are standard tools available
for any Linux distribution. Use you preferred package installer to have them available in your system. QoptKit also requires the C++ library Eigen 3 to be installed in your system.
This library is also available in most of the Linux distributions.</p>

<p align="justify">
Python libraries qiskit, matplotlib and numpy are required but they will be installed automatically. We recommend the use of conda but it is not required.
</p>

# 3 Documentation
The documentation can be found in https://soqcsadmin.github.io/QoptKIT/

# 4. Authorship #
<b>Javier Osca</b> <br>
soqcslib@gmail.com

<b>Jiri Vala</b> <br>
jiri.vala@mu.ie

# 5. License and Copyright #
Copyright (c) 2023 National University of Ireland Maynooth, Maynooth University. 

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the [License](./LICENSE.TXT). 

