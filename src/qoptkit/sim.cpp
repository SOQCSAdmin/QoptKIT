//======================================================================================================
// File sim.cpp
//
// SIMULATOR LIBRARY
//
// Copyright 2023 National University of Ireland Maynooth
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//=============================================================================

#include "sim.h"
#include <iostream>
#include <fstream>


//----------------------------------------
//
// Create a circuit "simulator".
//
//----------------------------------------
simulator::simulator(int i_mem){
//  int i_mem            // Number of memory positions reserved.


    mem=i_mem;
}


//----------------------------------------
//
// Destroy circuit "simulator"
//
//----------------------------------------
simulator::~simulator(){
}


//--------------------------------------------------------------
//
// Permanent calculation method (using Glynn formula). Full distribution.
//
//---------------------------------------------------------------
state *simulator::run( state *istate, qocircuit *qoc ){
//  state     *istate;           // Input state
//  qocircuit *qoc               // Circuit to be simulated
//  Variables
    int    nph;                  // Number of photons present in input ket.
    int    nlevel;               // Number of levels (qoc has this information, but it is put in this variable for easy access)
    int    index;                // Ket list position where a new term of the output state is stored.
    int   *pos;                  // Level where each photon is located. "Photon position"
    int   *occ;                  // Occupation
    cmplx  coef;                 // Coefficient for the transformation of a ket.
    cmplx  s;                    // Normalization coefficient of the input ket
    cmplx  t;                    // Normalization coefficient of the output ket
    matc   Ust;                  // Matrix to calculate the permanent
    state *ostate;               // Output state
//  Index
    int    iket;                 // Index of input kets elements
    int    ilin;                 // Index of input levels
    int    ilout;                // Index of output levels
    int    irow;                 // Row index of Ust
    int    icol;                 // Col index of Ust
//  Auxiliary index
    int    i;                    // Aux index
    int    j;                    // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(istate->nph,nlevel,mem);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate variables from the input ket
        nph=0;
        s=1.0;
        for(i=0;i<nlevel;i++){
            nph=nph+istate->ket[iket][i];
            s=s*(cmplx)factorial(istate->ket[iket][i]);
        }


        // Check all the possible outputs
        pos=new int[nph+1]();
        occ=new int[nlevel]();
        Ust.setZero(nph,nph);

        while(pos[0] < nlevel){
             // Calculate variables from the output ket
            t=1.0;
            for(j=0;j<nlevel;j++) occ[j]=0;
            for(j=0;j<nph;j++) {
                occ[pos[j]]=occ[pos[j]]+1;
                t=t*(cmplx)occ[pos[j]]; // This is the factorial implicitly.
            }


            // If the number of photons coincide (it always should)
            if(nph>0){
                // Create Ust
                icol=0;
                for(ilin=0;ilin<nlevel;ilin++){
                for(i=0;i<istate->ket[iket][ilin];i++){
                    irow=0;
                    for(ilout=0;ilout<nlevel;ilout++){
                    for(j=0;j<occ[ilout];j++){
                        Ust(irow,icol)=qoc->circmtx(ilout,ilin);
                        irow=irow+1;
                    }}
                    icol=icol+1;
                }}

                // Calculate coefficient
                coef=istate->ampl[iket]*glynn(Ust)/(sqrt(t)*sqrt(s));
            }else{
                coef=1.0;
             }


            // Store
            if(abs(coef)>xcut){
                index= ostate->add_term(coef,occ);
                if(index<0){
                    cout << "Simulator(GlynnF): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                    // Free memory
                    delete[] pos;
                    delete[] occ;
                    // Return partial calculation
                    return ostate;
                }
            }


            // Obtain new photon level "position"
            pos[nph-1] += 1; // xxxxN -> xxxxN+1
            for (i = nph; i > 0; i -= 1) {
                if (pos[i] > nlevel - 1) // if number spilled over: xx0(n-1)xx
                {
                    pos[i - 1] += 1; // set xx1(n-1)xx
                    for (j = i; j <= nph; j += 1)
                        pos[j] = pos[j - 1]; // set xx11..1
                }
            }
        }

        // Free memory
        delete[] pos;
        delete[] occ;

    }}
    // Return output
    return ostate;
}


//---------------------------------------------------------------------------
//
//  Assemble circuit. ( Path encoding version )
//  Creates a matrix relating the qubit encoded input and output states
//  to export to QISKIT.
//  WARNING: Unitarity is not guaranteed!!
//
//----------------------------------------------------------------------------
cmplx *simulator::assemble(mati in_qdef, mati out_qdef, int* occancilla, int *cond, qocircuit *qoc, int norm){
// mati        qdef;       // Integer matrix with the qubit to channel definitions
// state      *ancilla;    // We need an ancilla state to inform the decoding process of the auxiliary non-qubit channels values
// qocircuit  *qoc;        // Quantum optical circuit to which this state is referred
// int         method;     // Simulation method used to assemble the matrix
// int         norm;       // Is the output matrix normalized 0='No'/1='Yes'
// Variables
    int        nqb;        // Number of qubits
    int        ncoef;      // Number of coefficients
    int        maxnph;
    cmplx     *mtx;        // Assembled matrix. (In vector form to make communication with python easy)
    state     *qbits;      // Initial qubit encoded state
    state     *input;      // Initial photon encoded state
    state     *output;     // Photonic output state
    state     *pselected;  // Photonic output state after post-selection
    state     *encoded;    // Output qubit encoded state
    state     *ancilla;
    projector *prj;        // Projector
    long long  int ival;   // Matrix row index
    long long  int oval;   // Matrix col index
//  Auxiliary index
    int i;                 // Aux index
    int j;                 // Aux index


    // Initialize parameters
    nqb=in_qdef.cols();
    ncoef=pow(2,nqb);
    mtx=new cmplx[ncoef*ncoef]();

    // Check definition
    if(nqb!=out_qdef.cols()){
        cout << "Assemble error: Different number of qubits in the output with respect the input." << endl;
        return mtx;
    }

    //Initialize qubit input state
    qbits=new state(1,nqb+1,1);
    qbits->ampl[0]=1.0;
    qbits->nket=1;

    maxnph=in_qdef.cols();
    for(i=0;i<qoc->nlevel;i++) maxnph=maxnph+occancilla[i];
    ancilla=new state(maxnph,qoc->nlevel,1);
    ancilla->add_term(1.0,occancilla);

    // Initialize post-selection projector
    // We are limited to a single ket projector
    prj=new projector(maxnph,qoc->nlevel,1);
    prj->add_term(1.0,cond);

    // MAIN LOOP: For all possible input qubit configuration
    while(qbits->ket[0][nqb]==0){
        // Decode the input into a photons state
        input=qbits->decode(in_qdef,ancilla, qoc);

        // Run a simulation to obtain the output
        output=run(input,qoc);

        // Perform post-selection if needed
        pselected=output->post_selection(prj);

        // Encode an normalize the photonic output into qubits (path encoding)
        encoded=pselected->encode(out_qdef,qoc);
        if(norm==1) encoded->normalize();

        // Store each element at the proper matrix position
        ival=decval(qbits->ket[0],nqb,2);
        for(i=0;i<encoded->nket;i++){
            oval=decval(encoded->ket[i],nqb,2);
            mtx[ncoef*oval+ival]=encoded->ampl[i];
        }

        // Free memory
        delete input;
        delete output;
        delete pselected;
        delete encoded;

        // Calculate the next input state
        qbits->ket[0][0]=qbits->ket[0][0]+1;
        for(j=0;j<nqb;j++){
            if(qbits->ket[0][j]>=2){
                qbits->ket[0][j]=0;
                qbits->ket[0][j+1]=qbits->ket[0][j+1]+1;
            }
        }
    }

    // Free memory
    delete qbits;
    delete ancilla;
    delete prj;

    // Return matrix
    return mtx;
}
