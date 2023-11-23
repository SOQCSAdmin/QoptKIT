/**************************************************************************
* @file sim.h
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
* @title Simulator library
*
***************************************************************************/



#include "state.h"
#include <thread>
#include <future>


class simulator{
public:
    // Public variables
    int mem;                                   // Memory reserved for operations

    // Public functions
    // Management functions
    simulator(int i_mem);                      // Creates a circuit simulator.
    ~simulator();                              // Destroys a circuit simulator.

    // Simulation execution functions
    state *run(state *istate,qocircuit *qoc);  // Calculates the output state

   // Compilation methods
    cmplx *assemble(mati in_qdef, mati out_qdef, int* occanzilla, int *cond, qocircuit *qoc, int norm); // Assemble circuit matrix in qubit encoding
};
