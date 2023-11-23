//======================================================================================================
// File qocircuit.cpp
//
// OPTICAL CIRCUIT LIBRARY.
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

#include "qocircuit.h"


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch){
//  int i_nch;  // Number of channels


    nlevel=i_nch;
    if(nlevel<=0){
        cout << "Create circuit error #1: Number of channels has to be greater than zero." << endl;
        nlevel=1;
    }

    // Create and initialize circuit.
    // No elements in circuit is the identity matrix
    circmtx=matc::Identity(nlevel,nlevel);
}


//----------------------------------------
//
//  Destroy circuit
//
//----------------------------------------
qocircuit::~qocircuit(){
}


//----------------------------------------
//
//  Copy a circuit
//
//----------------------------------------
qocircuit *qocircuit::clone(){
//  Variables
    qocircuit* newcircuit;      // New circuit recipient of the copy


    //Reserve memory
    newcircuit= new qocircuit(nlevel);

    // Copy the circuit configuration
    newcircuit->circmtx=circmtx;

    //Return a pointer to the new circuit object.
    return newcircuit;
}


//---------------------------------------------------------
//
//  Returns the number of levels of this circuit
//
//---------------------------------------------------------
int qocircuit::num_levels(){
    return nlevel;
}


//----------------------------------------
//
//  Adds a custom gate to the circuit
//
//----------------------------------------
void qocircuit:: custom_gate(veci CH, matc U){
//  mati iodef;       // List of channels and polarization the define the input to the gate
//  matc U;           // Custom gate nxn matrix definitions
//  Variables
    int  nbmch;       // U square matrix row/column dimension
    int  ch1;         // Channel evaluated 1
    int  ch2;         // Channel evaluated 2
    matc oelement;    // Extended nlevelxnlevel custom gate definition
//  Auxiliary index.
    int  k;           // Aux index
    int  l;           // Aux index.


    // Definitions format adaptation
    nbmch=U.cols();

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=CH(k);

        for(l=0;l<nbmch;l++){
            ch2=CH(l);
            oelement(ch1,ch2)=U(k,l);
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;

}


//----------------------------------------
//
//  Adds a beamsplitter to the circuit
//
//----------------------------------------
void qocircuit:: beamsplitter(int i_ch1, int i_ch2, double d_theta, double d_phi){
//  int    i_ch1;             // Beamsplitter input channel 1
//  int    i_ch2;             // Beamsplitter input channel 2
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees
//  Constants
    const int nbmch=2;        // U square matrix row/column dimension
//  Variables
    double theta;             // Angle theta in radians
    double phi;               // Angle phi in radians
    matc   U(nbmch,nbmch);    // Beamsplitter 2x2 matrix definitions
    veci   V(nbmch);        // Channels to which U refers


    // Conversion to radians
    theta=d_theta*pi/180.0;
    phi  =d_phi*pi  /180.0;

    // 2x2 Matrix initialization
    U(0,0)= cos(theta);
    U(0,1)=-exp( jm*phi)*sin(theta);
    U(1,0)= exp(-jm*phi)*sin(theta);
    U(1,1)= cos(theta);

    // Channels
    V(0)  = i_ch1;
    V(1)  = i_ch2;

    // Compute gate
    return custom_gate(V,U);
}


//----------------------------------------
//
// General phase shifter with losses
//
//----------------------------------------
void qocircuit:: phase_shifter(int i_ch, double d_phi){
//  int   i_ch                // Phase shifter input channel
//  cmplx t                   // Phase shifter transmission amplitude of probability.
//  Constants
    const  int nbmch=1;       // U square matrix row/column dimension
//  Variables
    double phi;               // Aux index.
    matc   U(nbmch,nbmch);    // Waveplate 2x2 matrix definitions
    veci   V(nbmch);          // Channels to which U refers

    // Conversion to radians
    phi=d_phi*pi/180.0;

    // Matrix
    U(0,0)=exp(jm*phi);

    // Channel
    V(0)=i_ch;

    // Compute gate
    custom_gate(V,U);
}


//------------------------------------------------
//
// Rewires / swaps two channels (of the same mode)
//
//------------------------------------------------
void qocircuit:: rewire(int i_ch1,int i_ch2){
//  int    i_ch1;             // Channel 1
//  int    i_ch2;             // Channel 2
//  Constants
    const int nbmch=2;        // U square matrix row/column dimension
//  Variables
    matc   U(nbmch,nbmch);    // Beamsplitter 2x2 matrix definitions
    veci   V(1,nbmch);        // Channels to which U refers


    // 2x2 Matrix initialization
    U(0,0)= 0.0;
    U(0,1)= 1.0;
    U(1,0)= 1.0;
    U(1,1)= 0.0;

    // Channels
    V(0)  = i_ch1;
    V(1)  = i_ch2;

    // Compute gate
    custom_gate(V,U);
}


//------------------------------------------------
//
// Adds a gate using other circuit as the gate definition
//
//------------------------------------------------
int qocircuit:: add_gate(int *chlist, qocircuit *qoc){
//  veci chlist;     // List of channels where the gate is added
//  qocircuit *qoc;  // Quantum optical circuit
//  Variables
    matc   U;        // Gate matrix definition
    veci   W;        // Channels to which U refers
//  Auxiliary index
    int    i;        // Aux index


    // Copy gate definition
    U=qoc->circmtx;
    W.resize(qoc->nlevel);
    for(i=0; i<qoc->nlevel; i++) W(i)=chlist[i];

    // Create gate
    custom_gate(W,U);

    // Return success
    return 0;
}

