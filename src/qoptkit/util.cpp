//=============================================================================
// File util.cpp
//
// CONSTANTS and UTILITIES LIBRARY
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

#include "util.h"
#include <chrono>


// Initialization of extern variables
int def_nph=4;        // Default value of the maximum photon occupation by level.


//-----------------------------------------------
//
// Calculates the power p of an integer x.
// C++ pow only works with floats.
// Recipe obtained from:
// https://stackoverflow.com/questions/1505675/power-of-an-integer-in-c
// by Matthieu M. and edited by Johan
//
//-----------------------------------------------
int intpow(int x, unsigned int p)
{
//  int x;              // Base
//  unsigned int p;     // Power


  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = intpow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}


//--------------------------------------------
//
//  Calculate the factorial of n
//  Taken from literature/web. Various sources (common knowledge)
//
//--------------------------------------------
long int factorial(long int n){
//  long int n;      // Integer to calculate the factorial


    // Single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}


//--------------------------------------------
//
//  Sets the maximum number of photons by level in the
//  simulation and the debug flag. It establishes a base
//  for the indexing mechanism
//
//--------------------------------------------
void cfg_qoptkit(int nph){
//  int nph     // Number of photons by level


    def_nph=nph;
}


//--------------------------------------------
//
//  Calculate a hash value for a given vector of numbers in a
//  for a maximum number of photons. (We get an unique number for a given occupation vector)
//
//--------------------------------------------
long long int hashval(int *chainv,int n,int nph){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv
//  int nph           // Maximum number of photons


    return decval(chainv,n,nph+1);
}


//--------------------------------------------
//
//  Calculate a decimal value for a given vector of numbers in a
//  specified base.
//
//--------------------------------------------
long long int decval(int *chainv,int n,int base){
//  int *chainv       // List of numbers. For us occupations.
//  int n             // Length of chainv
//  int base          // Base of the numbers in chainv
//  Variables
    long long int value;   // Value to be returned
//  Auxiliary index
    int i;            // Aux index


    value=0;
    for(i=0;i<n;i++){
        if(chainv[i]>=0){
            value=value*base+chainv[i];
        }
    }
    return value;
}


//-----------------------------------------------
//
//  Calculation of the permanent of a matrix using
//  the Glynn method in Gray code. Translation to C++
//  of the algorithm in python coded by "xnor" in:
//  https://codegolf.stackexchange.com/questions/97060/calculate-the-permanent-as-quickly-as-possible
//
//  Glynn Formula:
//  https://en.wikipedia.org/wiki/Computing_the_permanent
//
//-----------------------------------------------
cmplx glynn(matc M){
//  matc M Square matrix to calculate the permanent.
//  Variables
//  Iteration variables
    long int n;         // Number of row/columns of the square matrix M.
    cmplx total;        // Total value of the permanent.
    long int old_gray;  // Old gray codification
    long int new_gray;  // New gray codification
    long int num_loops; // Number of loops
    long int diff;      // old_gray - new_gray
    long int gray_diff; // old_gray^new_grays
    long int bin_index; // Binary index
    cmplx sign;         // Sing value
    cmplx direction;    // Direction of the difference.
    cmplx reduce;       // Reduced row combination
    cmplx *row_comb;    // Row combination
    cmplx *new_vector;  // New vector
    long int gray_diff_index;                   // Gray code index of differences
    thash binary_power_dict;                    // Hash table of binary power differences.
    thash::const_iterator hash_gray_diff_index; // Iterator for the has table of binary power differences.
//  Auxiliary index.
    long int i;         // Aux index
    long int j;         // Aux index


    // Configuration
    n=M.cols();
    if(n==0) return 1.0;

    // Initializations
    row_comb=new cmplx[n];
    for(j=0;j<n;j++){
        row_comb[j]=0.0;
        for(i=0;i<n;i++){
            row_comb[j]=row_comb[j]+M(i,j);
        }
    }

    total = 0;
    old_gray = 0;
    sign = +1;

    new_vector=new cmplx[n];
    for(i=0;i<n;i++) binary_power_dict[pow(2,i)]=i;
    num_loops=pow(2,n-1);


    //  Main loop
    for(bin_index=1;bin_index<=num_loops;bin_index++){

        reduce=1.0;
        for(i=0;i<n;i++) reduce=reduce*row_comb[i];
        total=total+(sign*reduce);

        new_gray = bin_index^(bin_index/2);
        gray_diff = old_gray^new_gray;
        gray_diff_index = binary_power_dict[gray_diff];
        hash_gray_diff_index=binary_power_dict.find(gray_diff);
        gray_diff_index=hash_gray_diff_index->second;

        for(i=0;i<n;i++) new_vector[i]=M(gray_diff_index,i);
        diff=old_gray-new_gray;
        direction=0;
        if(diff>0) direction=2;
        if(diff<0) direction=-2;

        for(i=0;i<n;i++){
            row_comb[i] = row_comb[i]+ (new_vector[i] * direction);
        }

        sign = -sign;
        old_gray = new_gray;
    }

    // Free memory
    delete[] row_comb;
    delete[] new_vector;

    // Return value
    return total/(double)num_loops;

}
