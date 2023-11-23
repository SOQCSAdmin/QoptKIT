//======================================================================================================
// File state.cpp
//
// STATE DEFINITION LIBRARY.
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

#include "state.h"


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets.
//
//-----------------------------------------------------
ket_list::ket_list(int i_nph, int i_level, int i_maxket){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list


    create_ket_list(i_nph, i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
ket_list::ket_list(int i_nph, int i_level, int i_maxket, int *i_vis){
//  int  i_nph       // Maximum number of photons
//  int  i_level     // Number of levels to describe a ket
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels
//  Auxiliary  index
    int  g;           // Aux index


    create_ket_list(i_nph, i_level,i_maxket);
    for(g=0;g<i_level;g++) vis[g]=i_vis[g];
}


//----------------------------------------
//
//  Auxiliary private method to create a ket list
//  Each ket is distinguished by level occupations
//
//----------------------------------------
void ket_list::create_ket_list(int i_nph, int i_level, int i_maxket){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list
//  Auxiliary index
    int i;          // Aux index


    // Init ket information
    if(i_nph>0) nph=i_nph;
    else nph=def_nph;
    nket=0;
    nlevel=i_level;
    maxket=i_maxket;

    // Create and int amplitude/coefficient and ket structures
    ket=new int*[maxket];
    for(i=0;i<maxket;i++){
        ket[i]=new int[nlevel]();
    }

    // Compute the trivial print visibility vector.
    vis=new int[nlevel];
    for(i=0;i<nlevel;i++){
        vis[i]=i;
    }
}



//----------------------------------------
//
//  Destroys a ket_list
//
//----------------------------------------
ket_list::~ket_list(){
//  Auxiliary index
    int i;          // Aux index


    // Liberate memory of amplitude/coefficients and ket structures
    for(i=0;i<maxket;i++){
        delete[] ket[i];
    }

    //Free memory
    delete[] ket;
    delete[] vis;
    // Clear hash table
    ketindex.clear();
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
ket_list *ket_list::clone(){
    // Variable
    ket_list *aux;  // Auxiliary state recipient of the copy
    // Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    aux=new ket_list(nph, nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }
    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
// Finds the position of a ket in the list
//
//----------------------------------------
int ket_list::find_ket(int *occ){
//  int *occ;                    // Occupation of those levels
//  Variables
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.


    // Update amplitude and occupation
    // Store
    value=hashval(occ,nlevel,nph);
    hashv=ketindex.find(value);
    // If the ket does not exist yet on the output create new one.
    if(hashv==ketindex.end()){
        return -1;
    }else{
        return hashv->second;
    }
}


//----------------------------------------
//
//  Adds a new ket to the list state
//
//----------------------------------------
int ket_list::add_ket(int *occ){
//  int *occ;                    // Level occupation
//  Variables
    int  index;                  // Index/List position to be returned
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.
//  Auxiliary index
    int  i;                      // Aux index


    // Check and warn about memory limits.
    if(nket>=maxket){
        cout << "Ket list add ket error #1: Memory limit exceeded." << endl;
        return -1;
    }

    // Update amplitude and occupation
    //Store
    value=hashval(occ,nlevel,nph);
    hashv=ketindex.find(value);

    if(hashv==ketindex.end()){
        index=nket;
        ketindex[value]=nket;
        for(i=0;i<nlevel;i++){
            ket[nket][i]=occ[i];

        }

        nket=nket+1;
    }else{
        index=hashv->second;
    }

    // Return storage index
    return index;
}



//-------------------------------------------------------
//
//  Clear the list (without destroying it).
//
//-------------------------------------------------------
void ket_list::clear_kets(){


    nket=0;
    ketindex.clear();
}


//--------------------------------------------
//
//  Prints a ket
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket){
//  int iket;           // Ket number to be printed
//  int format;         // Integer that defines the notation used when printing the ket
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int writeprev;      // Have we print a level of this ket previously? 1=Yes/0=No
//  Auxiliary index
    int k;              // Aux index


    cout  << CYAN << " | ";
    writeprev=0;
    for(k=0;k<nlevel;k++){
        if(writeprev==1) cout << ", ";
        if(ket[iket][k]>=0){
            cout << ket[iket][k];
        }else{
            cout << "X";
        }
        writeprev=1;

    }
    cout << " >"<< RESET;
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets.
//
//-----------------------------------------------------
state::state(int i_nph, int i_level, int i_maxket):ket_list(i_nph, i_level,i_maxket){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe the state
//  int i_maxket    // Maximum number of different kets in the summation


    ampl=new cmplx[maxket]();
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
state::state(int i_nph, int i_level, int i_maxket, int *i_vis):ket_list(i_nph, i_level,i_maxket, i_vis){
//  int  i_nph       // Maximum number of photons
//  int  i_level     // Number of levels to describe the state
//  int  i_maxket    // Maximum number of different kets in the state
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    ampl=new cmplx[maxket]();
}


//----------------------------------------
//
//  Destroys a state
//
//----------------------------------------
state::~state(){


    // Liberate memory of amplitude
    delete[] ampl;
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
state *state::clone(){
    // Variable
    state *aux;  // Auxiliary state
    // Auxiliary index
    int    i;    // Aux index
    int    j;    // Aux index


    aux=new state(nph, nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        aux->ampl[i]=ampl[i];
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }
    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
//  Clear a state
//
//----------------------------------------
void state::clear(){


    delete[] ampl;
    ampl=new cmplx[maxket]();
    clear_kets();
}


//----------------------------------------
//
//  Adds a new term to a state
//
//----------------------------------------
int state::add_term(cmplx i_ampl, int *occ){
//  cmplx i_ampl;       // Amplitude of the new term
//  int  *occ;          // Occupation of those levels
//  Variables
    int   index;        // Index to be returned


    index=add_ket(occ);
    if(index>=0) ampl[index]=ampl[index]+i_ampl;

    return index;
}


//---------------------------------------------------
//
// Calculates the braket <prj|state>
//
//---------------------------------------------------
cmplx state::braket(state *bra){
//  state *bra;     // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    cmplx  rampl;   // Amplitude to be returned
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    rampl=0.0;
    // For each projector term
    for(i=0;i<bra->nket;i++){
            // Find the corresponding ket in the state
            j=find_ket(bra->ket[i]);
            if(j>=0){
                rampl=rampl+conj(bra->ampl[i])*ampl[j];
            }
    }

    // Return state.
    return rampl;
}


//---------------------------------------------------
//
// Normalizes the state
//
//---------------------------------------------------
void state::normalize(){
//  Variables
    double tot;     // Total/ Norm of the state
//  Auxiliary index
    int    i;       // Aux index


    tot=0.0;
    for(i=0;i<nket;i++){
        tot=tot+real(conj(ampl[i])*ampl[i]);
    }
    if(abs(tot)>xcut){
        for(i=0;i<nket;i++){
            ampl[i]=ampl[i]/sqrt(tot);
        }
    }

}


//--------------------------------------------
//
//  Post selection process of the output state
//
//--------------------------------------------
state *state::post_selection(state *prj){
//  state *prj;      // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    int    npost;    // Number of levels in which post selection is performed
    int    selected; // Is this state selected 0=No/1=Yes
    int   *islincl;  // Is level included? 0=No/1=Yes
    int   *occ;      // Occupation
    state *nstate;   // New post-selected state
// Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index


    // Determined the number of post-selected levels and which ones
    // will be included in nstate or not.
    npost=0;
    islincl=new int[nlevel];
    for(i=0;i<prj->nlevel;i++){
        if(prj->ket[0][i]>=0){
            islincl[i]=0;
            npost++;
        }else{
            islincl[i]=1;
        }
    }

    // Create new post-selected state / Reserve memory
    nstate=new state(nph, nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nstate->vis[k]=vis[l];
            k++;
        }
    }

    // For each projector term
    occ=new int[nlevel-npost]();
    for(i=0;i<prj->nket;i++){
        // Post-select each ket
        for(j=0;j<nket;j++){
            // Check selection condition
            selected=1;
            k=0;
            while((k<nlevel)&&(selected==1)){
                if ((ket[j][k]!=prj->ket[i][k])&&(prj->ket[i][k]>=0))  selected=0;
                k++;
            }

            // If is selected create the list of levels and
            // occupations for those not post-selected
            if(selected==1){
                k=0;
                for(l=0;l<nlevel;l++){
                    if(islincl[l]==1){
                        occ[k]=ket[j][l];
                        k++;
                    }
                }
                nstate->add_term(ampl[j]*conj(prj->ampl[i]),occ);
            }
        }
    }

    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return state.
    return nstate;
}


//----------------------------------------
//
// Returns the occupation of the indexed bin
// in string format.
//
//----------------------------------------
string state::tag(int index){
//  int index;             // Index of the bin
//  Variables.
    int len;               // Length of the auxiliary string
    int diff;              // Difference with the number of levels
    long long int value;   // Hash value of the ket in that idex.
    string aux;            // Auxiliary string. Conversion of value
    string temp;           // Temporary string=aux+0's at the left to match nlevel length
//  Auxiliary index
    int i;                 // Aux index


    // Perform conversion
    value=decval(ket[index],nlevel,10);
    aux=to_string(value);
    // Match the nlevel lenfth with 0's at the left
    len=aux.length();
    diff=nlevel-len;
    for(i=0;i<diff;i++) temp.append("0",1);
    temp.append(aux);
    // Return result
    return temp;
}


//--------------------------------------------
//
// General method to print a state in columns as
// a list ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(){
//  int format;         // Notation in which the state is printed
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int  firstline;     // Is this the first term? 1=Yes/0=No
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        if (abs(ampl[i])>DEFTHOLDPRNT){
            prnt_ket(i);


            firstline=0;
            cout << ": ";

            if(real(ampl[i])>=0){
                if(imag(ampl[i])>=0) cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << " j" << endl;
                else                 cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << " j" << endl;
            }else{
                if(imag(ampl[i])>=0) cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << " j" << endl;
                else                 cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << " j" << endl;
            }

        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//----------------------------------------------------------------------
//
//  Direct product of states defined in non coincident channels
//  (Not to be used in general therefore private)
//
//----------------------------------------------------------------------
int state::dproduct(state *rhs){
//  state *rhs      // State in the right hand side to do the d-product
//  Variables
    int    index;   // Ket list position where a new term of the dproduct.
    int   *occ;     // Occupation
    state *aux;     // Auxiliary state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index


    // Initialize
    aux=new state(this->nph, this->nlevel,this->maxket);
    aux=this->clone();
    this->clear();
    occ= new int[nlevel]();

    // Dproduct
    for(i=0;i<aux->nket;i++){
        for(j=0;j<rhs->nket;j++){
            for(k=0;k<aux->nlevel;k++){
                occ[k]=aux->ket[i][k]+rhs->ket[j][k];
            }
            index=add_term(aux->ampl[i]*rhs->ampl[j],occ);
            if(index<0){
                cout << "Dproduct: Warning! Memory exceeded. Operation cancelled." << endl;
                // Free memory
                delete[] occ;
                delete aux;
                return -1;
            }
        }
    }

    // Free memory
    delete[] occ;
    delete aux;
    return 0;
}


//---------------------------------------------------------------------------
//
//  Encode state into a qubit representation using path encoding
//
//----------------------------------------------------------------------------
state *state::encode(mati qdef,qocircuit *qoc){
//  mati       qdef; // Integer matrix with the qubit to channel definitions
//  qocircuit *qoc;  // Quatum optical circuit to which this state is referred
//  Variables
    bool   valid;    // Is the ket valid for codification true=Yes/false=No
    bool   print;    // A warning has been printed. True=Yes/False=No
    int    val0;     // value of channel 0
    int    val1;     // value of channel 1
    int    qval;     // qubit equivalent value of channels 0 and 1.
    int    idx;      // Index of the stored state.
    int    nvalid;   // Number of encoded kets
    int   *values;   // Vector of the values of each qubit
    state *qstate;   // qubit encoded state
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index
    int    m;        // Aux index
    int    n;        // Aux index


    // Reserve memory
    qstate=new state(1,qdef.cols(),maxket);


    // Encode each ket
    nvalid=0;
    print=false;
    for(i=0;i<nket;i++){
        values=new int[qdef.cols()]();
        valid=true;
        for(j=0;j<qdef.cols();j++){
            // Find the channels
            m=qdef(0,j);
            n=qdef(1,j);

            // Read the values
            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((vis[k]!=m)&&(k<nlevel)) k=k+1; // Not efficient but small search
            while((vis[l]!=n)&&(l<nlevel)) l=l+1; // Not efficient but small search
            val0=ket[i][k];
            val1=ket[i][l];

            // Encode the values
            qval=-1;
            if((val0==0)&&(val1==1)) qval=0;
            if((val0==1)&&(val1==0)) qval=1;
            if(qval<0) valid=false;
            values[j]=qval;
        }
        // If the resulting state is valid store it
        if(valid==true){
                idx=qstate->add_term(ampl[i],values);
                if((idx!=nvalid)&&(print==false)){
                    cout << "Encode: Warning! encoding leads to collision. Invalid result" << endl;
                    print=true;
                }
                nvalid=nvalid+1;
        }
        delete values;
    }

    // Return encoded state
    return qstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  path encoding. State version.
//
//----------------------------------------------------------------------------
state *state::decode(mati qdef,state *ancilla,qocircuit *qoc){
// mati       qdef; // Integer matrix with the qubit to channel definitions
// state *ancilla;  // We need an ancilla state to inform the decoding process of the auxiliary non-qubit channels values
// qocircuit *qoc;  // Quantum optical circuit to which this state is referred
// Variables
    int    val0;    // value of channel 0
    int    val1;    // value of channel 1
    int   *occ;     // Occupation vector
    state *phstate; // Photonic state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index
    int    l;       // Aux index
    int    m;       // Aux index
    int    n;       // Aux index


    // Reserve and initialize memory
    phstate=new state(ancilla->nph,ancilla->nlevel,maxket);

    for(j=0;j<phstate->nlevel;j++) phstate->vis[j]=ancilla->vis[j];

    // Decode each ket
    for(i=0;i<nket;i++){
        occ=new int[phstate->nlevel]();
        for(j=0;j<phstate->nlevel;j++) occ[j]=ancilla->ket[0][j];
        for(j=0;j<qdef.cols();j++){
            // Find the channels
            m=qdef(0,j);
            n=qdef(1,j);

            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((phstate->vis[k]!=m)&&(k<phstate->nlevel)) k=k+1; // Not efficient but small search
            while((phstate->vis[l]!=n)&&(l<phstate->nlevel)) l=l+1; // Not efficient but small search


            // Decode qubit values
            if(ket[i][j]==0){
                val0=0;
                val1=1;
            }else{
                val0=1;
                val1=0;
            }
            occ[k]=val0;
            occ[l]=val1;
        }
        // Store the decoded ket into a new state
        phstate->add_term(ampl[i],occ);
        delete occ;
    }

    // Return result
    return phstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  path encoding. Vector version.
//
//----------------------------------------------------------------------------
state *state::decode(mati qdef,veci ancilla,qocircuit *qoc){
//  mati       qdef;     // Integer matrix with the qubit to channel definitions
//  veci ancilla;        // Vector defining the values of the ancilla channels in order ( from smaller to larger channel number )
//  qocircuit *qoc;      // Quantum optical circuit to which this state is referred
//  Variables
    int maxnph;
    int   *isquch;       // Is a qubit channel? For every channel 0=Ancilla/1= Qubit channel
    state *anzstate;     // Ansatz state. State with the ansatz channel values defined.
    state *decoded;      // Decoded state.
    int *occ;
//  Auxiliary index
    int    i;            // Aux index
    int    k;            // Aux index



    maxnph=qdef.cols();
    for(i=0;i<ancilla.size();i++) maxnph=maxnph+ancilla(i);


    // Determine which channel is ansatz
    isquch=new int[qoc->nlevel]();
    occ=new int[qoc->nlevel]();

    for(i=0;i<qdef.cols();i++){
        isquch[qdef(0,i)]=1;
        isquch[qdef(1,i)]=1;
    }

    // Define the state
    k=0;
    for(i=0;i<qoc->nlevel;i++){
        if(isquch[i]==1){
            occ[i]=0;
        }else{
            // Is ancilla
            occ[i]=ancilla[k];
            k=k+1;
        }
    }

    // Create the ansatz state
    anzstate=new state(maxnph,qoc->nlevel,1);
    anzstate->add_term(1.0,occ);

    // Decode
    decoded=decode(qdef,anzstate,qoc);

    // Free memory.
    delete occ;
    delete isquch;
    delete anzstate;

    // Return state
    return decoded;
}


//-------------------------------------------------------
//
// Creates a projector specifying the maximum number of kets.
//
//-------------------------------------------------------
projector::projector(int i_nph, int i_level, int i_maxket):state(i_nph, i_level, i_maxket){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation


    create_projector(i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a projector specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
projector::projector(int i_nph, int i_level, int i_maxket, int *i_vis):state(i_nph, i_level, i_maxket, i_vis){
//  int i_nph       // Maximum number of photons
//  int i_level      // Number of levels to describe the projector
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    create_projector(i_level,i_maxket);
}


//-------------------------------------------------------
//
//  Auxiliary method to create a projector.
//  Same that to create a state but non-defined
//  occupation levels are not assumed to be 0 occupied.
//  They are initialized with a negative value instead.
//  This way these levels are later on ignored/not
// considered in the post-selection operation
//
//-------------------------------------------------------
void projector::create_projector(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation
//  Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    for(i=0;i<maxket;i++){
        for(j=0;j<nlevel;j++){
            ket[i][j]=-1;
        }
    }
}
