/**************************************************************************
* @file state.h
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
* @title Bosonic state library
*
***************************************************************************/


#include "qocircuit.h"


// Constant defaults
const double DEFTHOLDPRNT= 0.0001;    // Default amplitude magnitude threshold for printing.


class ket_list{
public:
    // Public variables
    int nph;               // Maximum number of photons
    int nket;              // Number of kets (C1|1>+C2|2>+...+Cn|nket>.
    int maxket;            // Maximum number of kets.
    int nlevel;            // Number of levels in each ket |0, 1, 2, ... nlevel>.

    // Ket list definition
    thash ketindex;        // Hash table of the dynamic dictionary of kets
    int **ket;             // Ket definitions. Level occupations of each ket/term
    int *vis;              // Correspondence vector. Position to level index.
                           // It stores to which level correspond each vector position.
                           // After post-selection it keeps track of the original level number.

    // Public methods
    // Management methods
    ket_list(int i_nph, int i_level,int i_maxket);               // Creates a ket list
    ket_list(int i_nph, int i_level, int i_maxket, int *i_vis);  // Creates a ket list
    ~ket_list();                                                 // Destroy a ket list
    ket_list *clone();                                           // Copy a ket list
    void clear_kets();                                           // Reset/empties ket list

    // State manipulation methods.
    int add_ket(int *occ);                                       // Adds a new ket to a ket list
    int find_ket(int *occ);                                      // Finds the position of a ket in the list ( given its occupation description )

    //Print methods
    void  prnt_ket(int iket);                                    // Prints a ket

protected:
    void create_ket_list(int i_nph, int i_level, int i_maxket);  // Create ket list auxiliary function

};

class state: public ket_list{
public:
    // Public variables
    complex<double> *ampl;                                       // Amplitudes of each ket/term

    // Public methods
    // Management methods
    state(int i_nph,int i_level,int i_maxket);                   //  Creates a state.
    state(int i_nph,int i_level, int i_maxket, int *i_vis);      //  Creates a state.
    ~state();                                                    //  Destroys a state
    state *clone();                                              //  Copy a state
    void clear();                                                //  Clear/empty a state


    // State manipulation methods
    int add_term(cmplx i_ampl, int *occ);                        // Adds a new term to a state
    int dproduct(state *rhs);                                    // Direct product of states defined in non coincident channels).
    cmplx braket(state *bra);                                    // Calculates the braket <bra|state>
    void normalize();                                            // Normalizes the state
    state *post_selection(state *prj);                           // Post select a state using a "projector"
    string tag(int index);                                       // Returns the "label" of a state
    void  prnt_state();                                          // Prints a state


    // Qubit codification methods.
    state *encode(mati qbits,qocircuit *qoc);                    // Encode from photonic to qubit representation
    state *decode(mati qdef,state *ancilla,qocircuit *qoc);      // Decode from qubit to photonic representation
    state *decode(mati qdef,veci ancilla,qocircuit *qoc);        // Decode from qubit to photonic representation
};


class projector : public state{
public:
    // Public methods
    // Management methods
    projector(int i_nph,int i_level,int i_maxket);               //  Creates a projector.
    projector(int i_nph,int i_level, int i_maxket, int *i_vis);  //  Creates a projector.

    // Auxiliary methods
    void create_projector(int i_level, int i_maxket);            // Create projector auxiliary function
};
