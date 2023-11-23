//======================================================================================================
// File pyqoptkit.cpp
//
// INTERFACE WITH PYTHON LIBRARY. C++ SIDE.
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


//----------------------------------------
//
//  Free memory of arrays
//
//----------------------------------------
void free_mem(char *mem){
    delete mem;
}


//----------------------------------------
//
//  Converter of array to matrix of integers
//
//----------------------------------------

mati to_matitr(int *int_array, int n, int m){
    int i;
    int j;
    mati cnv;

    cnv.resize(m,n);
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
           cnv(j,i)=int_array[i*m+j];
        }
    }
    return cnv;
}


//----------------------------------------
//
//  Interface with Python. C++/C side.
//
//----------------------------------------
extern "C" {
    //--------------------------------------------------------------------------------------------------------------------------
    // GENERAL
    // Interface support
    void free_ptr(char *mem){ free_mem(mem); }
    // Configure QOPTKIT
    void all_cfg_qooptkit(int nph){cfg_qoptkit(nph);}
    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // QOCIRCUIT
    // Management functions
    long int qoc_new_qocircuit(int i_nch){ return (long int)new qocircuit(i_nch);}
    void qoc_destroy_qocircuit(long int qoc){qocircuit* aux=(qocircuit*)qoc; delete aux; }
    int qoc_num_levels(long int qoc){qocircuit* aux=(qocircuit*)qoc; return aux->num_levels();}

    // Circuit elements
    //      Basic elements
    void qoc_beamsplitter(long int qoc, int i_ch1, int i_ch2, double theta, double phi){qocircuit* aux=(qocircuit*)qoc; aux->beamsplitter(i_ch1,i_ch2, theta, phi);}
    void qoc_phase_shifter(long int qoc, int i_ch, double phi){ qocircuit* aux=(qocircuit*)qoc; aux->phase_shifter(i_ch,phi);}
    void qoc_rewire(long int qoc,int i_ch1,int i_ch2){ qocircuit* aux=(qocircuit*)qoc; aux->rewire(i_ch1, i_ch2);}
    void qoc_add_gate(long int qoc1, int *chlist, long int qoc2){qocircuit *auxqoc1=(qocircuit *)qoc1;
                                                                        qocircuit *auxqoc2=(qocircuit *)qoc2;
                                                                        auxqoc1->add_gate(chlist,auxqoc2);
                                                                        }

    //--------------------------------------------------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------------------------------------------------
    // STATE
    // Management functions
    long int st_new_state(int i_nph, int i_level, int i_maxket){ return (long int)new state(i_nph,i_level,i_maxket);}
    void st_destroy_state(long int st){state* aux=(state*)st; delete aux; }

    // State manipulation methods.
    double *st_braket(long int st1,long int st2){ state *auxst1=(state*)st1;
                                                  state *auxst2=(state*)st2;
                                                  cmplx value=auxst1->braket(auxst2);
                                                  double *auxdouble=new double[2];
                                                  auxdouble[0]=real(value);
                                                  auxdouble[1]=imag(value);
                                                  return auxdouble;}

    void st_normalize(long int st){ state *auxst=(state*)st; auxst->normalize(); }
    void st_add_term(long int st, double rampl, double iampl, int *term){ state *auxst=(state*)st;
                                                                          cmplx ampl=rampl+jm*iampl;
                                                                          auxst->add_term(ampl,term); }


    long int st_post_selection(long int st, long int prj){ state* auxst=(state *)st; projector* auxprj=(projector*)prj; return (long int)auxst->post_selection(auxprj); }

    // Print methods
    int st_nkets(long int st,  int index){state *auxst=(state*)st; return auxst->nket;}
    int st_nlevels(long int st){state *auxst=(state*)st; return auxst->nlevel;}
    char *st_tag(long int st, int index){state *auxst=(state*)st; string value =auxst->tag(index); char* char_array=new char[value.length()+1]; strcpy(char_array, value.c_str()); return char_array;}
    double *st_campl(long int st,  int index){state *auxst=(state*)st;
                                              cmplx value=auxst->ampl[index];
                                              double *auxdouble=new double[2];
                                              auxdouble[0]=real(value);
                                              auxdouble[1]=imag(value);
                                              return auxdouble;}

    void  st_prnt_state(long int st){  state* auxst=(state*)st; auxst->prnt_state(); cout << flush;}

    // Qubit codification methods
    long int st_encode(long int st, int *qdef, int nqbits, long int qoc)     {state *auxst=(state *)st;
                                                                              qocircuit *auxqoc=(qocircuit*)qoc;
                                                                              mati imat=to_matitr(qdef,nqbits,2);
                                                                              state*auxencoded=auxst->encode(imat,auxqoc);
                                                                              return (long int) auxencoded;}

    long int st_decode(long int st, int *qdef, int nqbits, long int ancilla, long int qoc)     {state *auxst=(state *)st;
                                                                                                state *auxan=(state *)ancilla;
                                                                                                qocircuit *auxqoc=(qocircuit*)qoc;
                                                                                                mati imat=to_matitr(qdef,nqbits,2);
                                                                                                state*auxdecoded=auxst->decode(imat,auxan,auxqoc);
                                                                                                return (long int) auxdecoded;}


    // PROJECTOR
    long int prj_new_projector(int i_nph, int i_level, int i_maxket){ return (long int)new projector(i_nph, i_level,i_maxket);}
    void prj_destroy_projector(long int prj){projector* aux=(projector *)prj; delete aux; }

    // State manipulation methods.
    void prj_add_term(long int prj, double rampl, double iampl, int *term){ projector *auxprj=(projector *)prj;
                                                                                               cmplx ampl=rampl+jm*iampl;
                                                                                               auxprj->add_term(ampl,term);
                                                                                               }
    //--------------------------------------------------------------------------------------------------------------------------



    //--------------------------------------------------------------------------------------------------------------------------
    // SIMULATOR
    // Management methods
    long int sim_new_simulator(int i_mem){ return (long int) new simulator(i_mem);}
    void sim_destroy_simulator(long int sim){simulator *aux=(simulator *)sim; delete aux; }

    // Run methods
    long int sim_run(long int sim,long int st,long int qoc ){ simulator *auxsim=(simulator *) sim; state  *auxst=(state *) st; qocircuit *auxqoc=(qocircuit*)qoc; return (long int) auxsim->run(auxst,auxqoc);}
    //--------------------------------------------------------------------------------------------------------------------------

   //--------------------------------------------------------------------------------------------------------------------------
   // QISKIT



    double *qk_assemble(long int sim, int *in_qdef, int *out_qdef, int nqbits, int *occancilla, int *cond, long int qoc, int norm) {
                                                                             simulator *auxsim=(simulator *)sim;
                                                                             qocircuit *auxqoc=(qocircuit*)qoc;
                                                                             mati imat=to_matitr(in_qdef,nqbits,2);
                                                                             mati omat=to_matitr(out_qdef,nqbits,2);
                                                                             cmplx *auxmtx=auxsim->assemble(imat,omat,occancilla,cond,auxqoc,norm);
                                                                             return (double *) auxmtx;}




   //--------------------------------------------------------------------------------------------------------------------------
}
