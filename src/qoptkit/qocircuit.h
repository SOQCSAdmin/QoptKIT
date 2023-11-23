/***********************************************************************************
* @file qocircuit.h
* @author Javier Osca
* @author Jiri Vala
*
* @copyright Copyright 2023 National University of Ireland Maynooth
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
* @title Optical circuit library
*
***********************************************************************************/


#include "util.h"


class qocircuit{
public:
    // Public variables
    int    nlevel;          // Number of levels.
    matc   circmtx;         // Matrix of the circuit as a function of the level number.


    // Public functions
    // Management functions
    qocircuit(int i_nch);   // Create circuit
    ~qocircuit();           // Destroy circuit
    qocircuit *clone();     // Copy a circuit
    void reset();           // Reset circuit
    int num_levels();       // Returns the number of levels of this circuit

    // Circuit elements
    void custom_gate(veci CH, matc U);         // Adds a custom gate
    void beamsplitter(int i_ch1, int i_ch2, double theta, double phi); // Adds a beamsplitter to the circuit
    void phase_shifter(int i_ch, double phi);  // Adds a phase shifter to the circuit
    void rewire(int i_ch1,int i_ch2);          // Adds a swap gate between two channels
    int add_gate(int *chlist, qocircuit *qoc); // Adds a new gate defined by a circuit
};
