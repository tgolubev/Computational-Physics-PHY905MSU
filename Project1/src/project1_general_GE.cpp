//Partb of Program1. This performs a general Gaussian Elimination.
//Imputs: The program requires 2 imputs from the command line: output_file_name n
//I.e. if we pass the program: out 3  it will repeat the algorithm for 10, 100, and 1000 mesh points.
//Computation time is calculated for each value of n (10^n = # of mesh points) seperately.

//Coded by: Timofey Golubev, Hao Lin, and Xingze Mao

//Include statements for header files
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>    //  Mathematical functions like sin, cos etc
#include <string>   //  useful library to operate on characters
#include <chrono>  //for high resolution clock

// use namespace std for output and input so don't need to write std:: in front of commands like cout, cin
using namespace std;
using namespace std::chrono; //for high res clock

// object for output files, for input files use ifstream
ofstream ofile;  //allows to work with output files in compact way

// Define functions used ahead of the main program,
//Inline function: the compiler takes the return statements 'return  1-.0*exp(-10.0*x)' and pastes them in the corect spots
//instead of calling the function every time.  This saves flops. Used for short (i.e. 1 line funcitons).
inline double ff(double x){return 100.0*exp(-10.0*x); //input is x in format double, returns: 100.0*exp(-10.0*x)
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program: returns 0 upon sucess, 1 if something went wrong --> define main program as an 'int' to do this. If don't want to do this
//can skip the 'int', can just write main(...
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
      double adiag_ratio;
      // Use 'new' for dynamic memory allocation. For n points we have n+1 elements (arrays start from 0)
      double *diagonal = new double[n+1];  //array  for diagonal matrix elements
      double *a = new double[n+1];  //array for non-diagonal matrix elements, not needed for the special case
      double *c = new double[n+1];  //array for non-diagonal elements
      double *f = new double[n+1];  //array for right hand side (known fnc f(x))
      double *x = new double[n+1];  //array for x-values
      double *v = new double[n+1];  //array for solution

      //start clock timer
      high_resolution_clock::time_point start = high_resolution_clock::now();

      //Dirichlet boundary  conditions
      v[0] = v[n] = 0.0;

      //Define x and f(x) values
      for (int i = 0; i <=n; i++){
          x[i] = i*h;
          f[i] = hh*ff(i*h);
      }

      //Define values of diagonal and non-diagonal matrix elements
      for (int i = 1; i<n; i++) {
          diagonal[i] = 2.0;
          a[i] = -1;
          c[i] = -1;
      }

      //Forward substition
      for (int i = 2; i <n; i++) {
          adiag_ratio = a[i-1]/diagonal[i-1];
          diagonal[i] = diagonal[i] - adiag_ratio*c[i-1];
          f[i] = f[i] - adiag_ratio*f[i-1];
      }

      //Backward substitution
      v[n-1] = f[n-1]/diagonal[n-1]; //1st linear eqn
      for (int i = n-2; i >0; i--) v[i] = (f[i]-v[i+1]*c[i])/diagonal[i];

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
         double RelativeError = fabs((exact(xval)-v[i])/exact(xval)); //fabs is absolute value
         ofile << setw(15) << setprecision(8) << xval;  //setprecision(8) sets that use 8 sigfigs
         ofile << setw(15) << setprecision(8) << v[i];
         ofile << setw(15) << setprecision(8) << exact(xval); //is analytical solution defined at top
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] diagonal; delete [] f; delete [] v; delete[] a; delete[] c; //memory that was occupied by these arrays is now freed

    }

    return 0;   // returns zero upon success, this is not required but can be useful.
}


