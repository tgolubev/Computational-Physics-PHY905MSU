//Project2 part d Version 1: The interacting case.
//Performs Jacobi rotation to find the eigenvalues of matrix which is formed from the discretization of the Schrodinger equation
//for 2 interacting electrons in a harmonic oscillator potential. The Coulomb interaction
//and potential they are in is described by vector V. Calculations are performed for different values of oscillator 'frequencies' w_r
//(describes strenth of the potential) specified in vector w_r_vector.
//The eigenvalues are also calculated by armadillo's eig_sym function to ensure accuracy of the Jacobi algorithm

//Input: The program requires the desired name of the output file and the number of mesh points to use for discretization.
//(i.e. meshpts=200-300 works well for the simple harmonic oscillator potential). An example input would be "out 200".

//Output: The groundstate eigenvectors are output to a single output file as a function of relative coordinate r.

//Coded by: Tim Golubev, Hao Lin, Xingze Mao

//Include statements for header files
#include <iostream>
#include <fstream>
#include <iomanip>  //for output file manipulation
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>    //  Mathematical functions like sin, cos etc
#include <string>   //  useful library to operate on characters
#include <chrono>  //for high resolution clock
#include <armadillo.h>

using namespace std;
using namespace arma;
using namespace std::chrono; //for high res clock

// object for output files, for input files use ifstream
ofstream ofile;  //allows to work with output files in compact way
ofstream eigval;

double offdiag(mat A, int& p, int& q, int n);   //pass by ref p and q
void Jacobi_rotate ( mat& A, mat& R, int k, int l, int n );

int main(int argc, char *argv[])
{

  int meshpts;
  string filename;
  string eigval_filename;
  //If in command-line has less than 1 argument, write out an error.
  if( argc <= 3 ){
       cout << "Bad Usage: " <<
            " read also file name on same line and number of mesh points" << endl;
        exit(1);
   }
   else{
      filename = argv[1];
      eigval_filename = argv[2];
      meshpts = atoi(argv[3]);  //atoi: convert ascii input to integer. Input number of mesh points to use.
   }

  // Declarations
  double rho_max = 50.0;  //problem will be solved over range [0, rho_max]. Range should be changed for using different w_r's.
  double alpha = 1; //set alpha=1
  vec w_r_vector;
  w_r_vector << 0.01 << 0.5 << 1.0 << 5 << endr;        //Define the w_r values we want to study
  int num_eig = 3;      //num of min eigvalues want to output

  double h = rho_max/((double) meshpts);  //define h=step size
  double hh = h*h;
  int n = meshpts -1;    //we skip the endpts
  vec rho(n);            //rho values vector ("x"-values)
  vec V(n);              //vector for the potential

  mat Output(n, 5);                 //Matrix for output of the eigenvectors for 4 different w_r's
  mat Output_eigvalues(num_eig,4);  //Matrix for output of the lowest eigenvalues for the 4 different w_r's

  //Fill vector rho (it stays the same for all w_r)
  for(int i=0; i<n; i++){
      rho(i) = (i+1)*h;         //1st element of rho will = h. We exclude the endpoints rho(0)=0 and rho(rho_max).
  }
  Output.col(0) = rho*alpha;  //fill 1st column with r=rho*alpha values

  //Main loop to repeat the eigenvalue computation for different values of w_r
  for(int k = 0; k<4; k++){
       double w_r = w_r_vector(k);                  //for each iteration, pick out one of the w_r values from the w_r vector
       cout << "w_r = " << w_r <<endl;
      //Fill potential vector V
      for(int i=0; i<n; i++){
          V(i) = rho(i)*rho(i)*w_r*w_r + 1/rho(i);  //harmonic oscillator potential with repulsive Coulomb interaction btw. 2 electrons
          //V(i) = rho(i)*rho(i)*w_r*w_r;           //harmonic oscillator potential WITHOUT repulsive Coulomb interaction (for testing)
       }

      //Define matrix
      mat A(n,n);
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
          if (i==j)
          {A(i,j) = 2.0/hh + V(i);
          }else if (abs(i-j)==1)
          {A(i,j) = -1.0/hh;
          }else {
              A(i,j) = 0.0;
          }
        }
      }

      //-----------------------------------------------------------------------------------------------
      // Testing with eig_sym  solver
      //start clock timer
      high_resolution_clock::time_point start1 = high_resolution_clock::now();
      vec eigval = eig_sym(A);
      //stop clock timer and output time duration
      high_resolution_clock::time_point finish1 = high_resolution_clock::now();
      duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
      cout << "Armadillo eig_sym CPU time = " << time1.count() << endl;
      //Output the lowest eigvalues from eig_sym solver
      vec lowest_eigval(num_eig);
      for (int i=0;i<num_eig;i++){
          lowest_eigval(i) = eigval(i);
      }
      cout<< lowest_eigval<<endl;
      //-----------------------------------------------------------------------------------------------

      //Jacobi rotation algorithm

      mat R = eye(n,n);       //initialize matrix for eigenvectors as the identity matrix
      double tolerance = 1.0E-10;
      int iterations = 0;
      double maxnondiag = 1.0;    //initialize
      int maxiter = n*n;          //anything <n^2 gives results that don't match eig_sym solver

      //start clock timer
      high_resolution_clock::time_point start2 = high_resolution_clock::now();
      //Perform Jacobi rotations always acting on the largest (abs value) non-diag element
      while ( maxnondiag > tolerance && iterations <= maxiter)
      {
         int p, q;
         maxnondiag  = offdiag(A, p, q, n);
         Jacobi_rotate(A, R, p, q, n);
         iterations++;
      }
      //stop clock timer and output time duration
      high_resolution_clock::time_point finish2 = high_resolution_clock::now();
      duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
      cout << "Jacobi rotate CPU time = " << time2.count() << endl;
      //cout<< maxnondiag << endl;

      //Vector with eigenvalues read off of diagonal elements of A
      vec eig_values(n);
      for (int i=0; i<n; i++){
         eig_values(i) = A(i,i);
      }

      //Find and output lowest eigenvalues
      double max_value = max(eig_values);
      int index;
      int groundstate_index;
      vec min_eig_values(num_eig);                                //vector to store min eigvalues
      vec groundstate_eigvector;
      for (int i=0; i<num_eig; i++){
           min_eig_values(i) = min(eig_values);
           index = index_min(eig_values);                         //find index of minimum eig_value
           if (i==0){
               groundstate_index = index;
               groundstate_eigvector = R.col(groundstate_index);  //arma synthax: extract column vector
           }
           eig_values(index) = max_value;  //set the minimum value found to = max value in matrix
       }
       cout<<min_eig_values<<endl;
       //Fill in output matrix with eigvector for current w_r
       Output.col(k+1) = groundstate_eigvector;
       //Fill in output eigvalues matrix
       Output_eigvalues.col(k) = min_eig_values;
  }

  //Setup output file
  ofile.open(filename);
  ofile << setiosflags(ios::showpoint | ios::uppercase); //sets to write i.e. 10^6 as E6

  // Write to file:
  double bndry_value = 0.0;

  ofile << "   w_r =        " << setprecision(4)<< w_r_vector(0);
  ofile<<"      "<< setprecision(4)<<w_r_vector(1) <<"       "<<setprecision(4)<<w_r_vector(2)<<"        "<<setprecision(4)<<w_r_vector(3)<<endl;
  ofile << "Rel. Coord.(r) " << "          Ground State Wavefunctions" <<endl;
  ofile << setprecision(8) << "   " <<bndry_value << "   " << bndry_value << "   "<<bndry_value <<"   " <<bndry_value <<"   "<< bndry_value<< endl;
  ofile << setprecision(8) << Output;
  ofile << setprecision(8) << "   " <<rho_max << "    " << bndry_value << "   "<<bndry_value <<"   " <<bndry_value <<"   "<< bndry_value<< endl;
  ofile.close();

  eigval.open(eigval_filename);
  eigval << setiosflags(ios::showpoint | ios::uppercase); //sets to write i.e. 10^6 as E6
  eigval << "   w_r =        " << endl << setprecision(4)<< "    " << w_r_vector(0);
  eigval <<"   "<< setprecision(4)<<w_r_vector(1) <<"    "<<setprecision(4)<<w_r_vector(2)<<"     "<<setprecision(4)<<w_r_vector(3)<<endl;
  eigval <<"                Eigenvalues" << endl;
  eigval << Output_eigvalues;
  eigval.close();

  return 0;


 }

//--------------------------------------------------------------------------------------------------------------------------------
//Functions

//  the offdiag function, using Armadillo

double offdiag(mat A, int& p, int& q, int n)       //pass by referenc p and q, allows them to be changed by fnc.
{
double max = 0.0;
for (int i = 0; i < n; ++i)
{
 for ( int j = i+1; j < n; ++j)
 {
     double aij = fabs(A(i,j));
     if ( aij > max)
     {
        max = aij;  p = i; q = j;
     }
 }
}
return max;
}


void Jacobi_rotate ( mat& A, mat& R, int k, int l, int n )
{
double s, c;
if ( A(k,l) != 0.0 ) {
double t, tau;
tau = (A(l,l) - A(k,k))/(2*A(k,l));

if ( tau >= 0 ) {
t = 1.0/(tau + sqrt(1.0 + tau*tau));
} else {
t = -1.0/(-tau +sqrt(1.0 + tau*tau));
}

c = 1/sqrt(1+t*t);
s = c*t;
} else {
c = 1.0;
s = 0.0;
}
double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
a_kk = A(k,k);
a_ll = A(l,l);
A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
A(l,k) = 0.0;  // same here
for ( int i = 0; i < n; i++ ) {
if ( i != k && i != l ) {
a_ik = A(i,k);
a_il = A(i,l);
A(i,k) = c*a_ik - s*a_il;
A(k,i) = A(i,k);
A(i,l) = c*a_il + s*a_ik;
A(l,i) = A(i,l);
}
//  And finally the new eigenvectors
r_ik = R(i,k);
r_il = R(i,l);

R(i,k) = c*r_ik - s*r_il;
R(i,l) = c*r_il + s*r_ik;
}
return;
} // end of function jacobi_rotate

