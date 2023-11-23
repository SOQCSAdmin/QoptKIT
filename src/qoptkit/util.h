/**************************************************************************
* @file util.h
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
* @title Constants and utilities library.
*
***************************************************************************/

//

//=============================================================================


#include <iostream>      // In/out
#include <iomanip>       // Formatting in/out
#include <complex>       // Complex numbers
#include <string>        // Strings management
#include <algorithm>     // Permutations
#include <unordered_map> // Hash tables
#include <eigen3/Eigen/Dense>   // Eigen3 library

// Name spaces
using namespace std;
using namespace Eigen;

//Types definitions.
typedef complex<double>  cmplx;            // Complex double precision definition.
typedef MatrixXi  mati;                    // Simplified type name of integer matrix.
typedef MatrixXd  matd;                    // Simplified type name of double matrix.
typedef MatrixXcd matc;                    // Simplified type name of complex double matrix.
typedef VectorXi  veci;                    // Simplified type name of integer vector.
typedef VectorXd  vecd;                    // Simplified type name of double vector.
typedef VectorXcd vecc;                    // Simplified type name of complex double vector.

// Hash table definition
typedef unordered_map<long long int,long long int> thash; // Simplified type name for a hash table


//Extern variables.
extern int def_nph;                        // Default value of the maximum photon occupation by level.

// Constant to be used across the library
const double xcut = 1.0e-10;               // Value below which a real number is truncated to zero.
const double pi   = std::acos(-1);         // Value of Pi.
const std::complex<double> jm(0, 1);       // The pure imaginary number i.


// Color codes for input in bash terminal.
// Taken from literature web (common knowledge).
// the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"                  // RESET.
#define BLACK   "\033[30m"                 // Black.
#define RED     "\033[31m"                 // Red.
#define GREEN   "\033[32m"                 // Green.
#define YELLOW  "\033[33m"                 // Yellow.
#define BLUE    "\033[34m"                 // Blue.
#define MAGENTA "\033[35m"                 // Magenta.
#define CYAN    "\033[36m"                 // Cyan.
#define WHITE   "\033[37m"                 // White.
#define BOLDBLACK   "\033[1m\033[30m"      // Bold Black.
#define BOLDRED     "\033[1m\033[31m"      // Bold Red.
#define BOLDGREEN   "\033[1m\033[32m"      // Bold Green.
#define BOLDYELLOW  "\033[1m\033[33m"      // Bold Yellow.
#define BOLDBLUE    "\033[1m\033[34m"      // Bold Blue.
#define BOLDMAGENTA "\033[1m\033[35m"      // Bold Magenta.
#define BOLDCYAN    "\033[1m\033[36m"      // Bold Cyan.
#define BOLDWHITE   "\033[1m\033[37m"      // Bold White.


//**************************************************************************
//
// Extra functions
//
//*************************************************************************

int intpow(int x, unsigned int p);                // Power p of an integer x.
long int factorial(long int n);                   // Factorial of integer n.
void cfg_qoptkit(int nph);                          // Configure the default number of photons.
long long int hashval(int *chainv,int n,int nph); // Hash value calculation of an occupation vector
long long int decval(int *chainv,int n,int base); // Decimal value for a given vector treating each entry as a digit.
cmplx glynn(matc M);                              // Permanent of a square matrix using the Balasubramanian/Bax/Franklin/Glynn formula in gray code (this is a compact and fast bit representation).
