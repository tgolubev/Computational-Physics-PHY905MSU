//LU decomposition for Program1.
//Imputs: The program requires 2 imputs from the command line: output_file_name n
//I.e. if we pass the program: out 3  it will repeat the algorithm for 10, 100, and 1000 mesh points.
//Computation time is calculated for each value of n (10^n = # of mesh points) seperately.

//Coded by: Timofey Golubev, Hao Lin, and Xingze Mao

//Include statements for header files
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>    //  Mathematical functions like sin, cos etc
#include <string>   //  useful library to operate on characters
#include <chrono>  //for high resolution clock

// use namespace std for output and input so don't need to write std:: in front of commands like cout, cin
using namespace std;
using namespace std::chrono; //for high res clock
#define   ZERO       1.0E-15

double ** AllocateMatrix(int, int);     //** double pointers are used for matrices
void DeallocateMatrix(double **, int, int);
//void WriteMatrix(double **, int);
//void MatrixMultiplication(double **, double **, int);
void LUDecomposition(double **, int, int *);
void LUBackwardSubstitution(double **, int, int *, double *);

// object for output files, for input files use ifstream
ofstream ofile;  //allows to work with output files in compact way

// Define functions used ahead of the main program,
//Inline function: the compiler takes the return statements 'return  1-.0*exp(-10.0*x)' and pastes them in the corect spots
//instead of calling the function every time.  This saves flops. Used for short (i.e. 1 line funcitons).
inline double ff(double x){return 100.0*exp(-10.0*x); //input is x in format double, returns: 100.0*exp(-10.0*x)
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program:
int main(int argc, char *argv[]){
  int exponent; 
  string filename; //  This is a useful way to define a text string, alternatively you can use characters

    //If in command-line has less than 1 argument, write out an error.
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1];        // first command line argument after name of program is put into element 1 of array: we use this for name of output file
        exponent = atoi(argv[2]);  //atoi: convert ascii input to integer. This is max power of 10^n.
    }

    // Loop over powers of 10 to allow to rerun the calculations for different n values, n = # of mesh (aka integration) pts
    //So this way can right away perform same calculation with different # of mesh pts. i.e. if i = 1, n = 10.  if i = 3, n=1000 etc.
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);         //pow(10.0) does 10^i and pow returns float by default. use pow(10.0,i) to return an int.
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-:  //when run the program (i.e. if typed 'test.exe out 10', it will make filenames be out1, out2 etc..
      fileout.append(argument);          //appends the 'argument' to our output file filename. The arguments will be i (i.e. 1, 2, 3)

      //   Declare Variables, Pointers, and Arrays
      double h = 1.0/((double) n);  //define h=step size. 1.0 to make sure it returns a double. write: ((double) n) to make sure n is double also.
      double hh = h*h;              //define h^2 as hh
      // Use 'new' for dynamic memory allocation. For n points we have n+1 elements (arrays start from 0)
      double *f = new double[n];  //array for right hand side (known fnc f(x))
      double *x = new double[n];  //array for x-values
      //For LU specifically
      double **A = AllocateMatrix(n, n);
      int *indx = new int[n];

      //Dirichlet boundary  conditions
      //f[0] = f[n] = 0.0;

      //start clock timer
      high_resolution_clock::time_point start = high_resolution_clock::now();

      //Define x and f(x) values
      for (int i = 0; i <n; i++){
          x[i] = i*h;
          f[i] = hh*ff(i*h);
      }

      //Define matrix for LU decomp
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
          if (i==j)
          {A[i][j] = 2.0;
          }else if (abs(i-j)==1)
          {A[i][j] = -1.0;
          }else {
              A[i][j] = 0.0;
          }
        }
      }

      // Perform the LU decomposition and backward substitution. Returns soln. in array f.
      LUDecomposition(A, n, indx);   // LU decompose  a[][]
      LUBackwardSubstitution(A, n, indx, f);

      //stop clock timer and output time duration
      high_resolution_clock::time_point finish = high_resolution_clock::now();
      duration<double> time = duration_cast<duration<double>>(finish-start);
      cout << "time = " << time.count() << endl;

      //Setup output file
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase); //sets to write i.e. 10^6 as E6

      //Write to file:
      //     ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
         double RelativeError = fabs((exact(xval)-f[i])/exact(xval)); //fabs is absolute value
         ofile << setw(15) << setprecision(8) << xval;   //setprecision(8) sets that use 8 sigfigs
         ofile << setw(15) << setprecision(8) << f[i-1]; // to correctly reshift index
         ofile << setw(15) << setprecision(8) << exact(xval); //is analytical solution defined at top
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] f; delete []indx; //memory that was occupied by these arrays is now freed
      DeallocateMatrix(A, n, n);     // release local memory
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------
//Matrix manipulation functions

// Allocate memory for a matrix and initialize the elements to zero
double ** AllocateMatrix(int m, int n){
  double ** Matrix;
  Matrix = new double*[m];
  for(int i=0;i<m;i++){
    Matrix[i] = new double[n];
    for(int j=0;j<m;j++)
      Matrix[i][j] = 0.0;
  }
  return Matrix;
}

// Write out a given matrix
void WriteMatrix(double ** Matrix, int n){
  for(int i=0;i < n;i++){
    cout << endl;
     for (int j=0 ; j < n;j++){
        printf("  A[%2d][%2d] = %12.4E",i, j, Matrix[i][j]);
     }
  }
    cout << endl;
}

// Free memory
void DeallocateMatrix(double ** Matrix, int m, int n){
  for(int i=0;i<m;i++)
    delete[] Matrix[i];
  delete[] Matrix;
}

// Straightforward matrix-matrix multiplication
/*
void MatrixMultiplication(double ** a, double **b, int n){
  double **c = AllocateMatrix(n, n);
  for(int i=0;i < n;i++){
     for (int j=0 ; j < n;j++){
       double sum = 0.0;
       for (int k = 0; k < n; k++) sum += a[i][k]*b[k][j];
       c[i][j] = sum;
     }
  }
  WriteMatrix(c,n);
}
*/

/*
    The function
    void LUDecomposition(double **a, int n, int *indx)
    takes as input a two-dimensional matrix a[][] of dimension n and
    replaces it by the LU decomposition of a rowwise permutation of
    itself. The results is stored in a[][]
    The vector
    indx[] records the row permutation effected by the partial pivoting.
*/

void LUDecomposition(double **a, int n, int *indx)
{
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new double [n];
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a[i][j])) > big) big = temp;
      }
      if(big == ZERO) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }
      vv[i] = 1.0/big;
   }

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a[i][j];
         for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a[i][j];
         for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         if((dum = vv[i]*fabs(sum)) >= big) {
            big = dum;
            imax = i;
         }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
            dum        = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k]    = dum;
         }
         vv[imax] = vv[j];         // also interchange scaling factor
      }
      indx[j] = imax;
      if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;
      if(j < (n - 1)) {                   // divide by pivot element
         dum = 1.0/a[j][j];
         for(i=j+1;i < n; i++) a[i][j] *= dum;
      }
   } // end j-loop over columns
   delete [] vv;   // release local memory
}


/*
     The function
       void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
     solves the set of linear equations A X = B of dimension n.
     a[][] is input, not as the matrix A[][] but rather as
     its LU decomposition, indx[] is input as the permutation vector returned by
     ludcmp(). b[] is input as the right-hand side vector B,
     The solution X is returned in B. The input data a[][],
     n and indx[] are not modified. This routine takes into
     account the possibility that b[] will begin with many
     zero elements, so it is efficient for use in matrix
     inversion.
*/

void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
{
   int        i, ii = -1, ip, j;
   double     sum;

   for(i = 0; i< n; i++) {
      ip    = indx[i];
      sum   = b[ip];
      b[ip] = b[i];
      if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
      else if(sum) ii = i;
      b[i] = sum;
   }
   for(i = n - 1; i >= 0; i--) {
      sum = b[i];
      for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
      b[i] = sum/a[i][i];
   }
}


