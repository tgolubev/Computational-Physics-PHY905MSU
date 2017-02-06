//Partc of Program1
//Include statements for header files
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>    //  Mathematical functions like sin, cos etc
#include <string>   //  useful library to operate on characters
#include <chrono>  //for high resolution clock

// use namespace for output and input so don't need to write std:: in front of commands like cout, cin
using namespace std;
using namespace std::chrono; //for high res clock

// object for output files, for input files use ifstream
ofstream ofile;  //allows to work with output files in compact way

// Define functions used ahead of the main program.
//Inline function: the compiler takes the return statements 'return  1-.0*exp(-10.0*x)' and pasts them in the corect spots
//instead of calling the function every time. This saves flops. Used for short (i.e. 1 line funcitons).

inline double ff(double x){return 100.0*exp(-10.0*x); //input is x in format double. returns: 100.0*exp(-10.0*x)
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program: returns 0 upon sucess, 1 if something went wrong --> define main program as an 'int' to do this. If don't want to do this
//can skip the 'int', can just write main(...
int main(int argc, char *argv[]){
  int exponent; 
  string filename; //  This is a useful way to define a text string, alternatively you can use characters

  //If in command-line has less than 1 argument, it writes out an error.
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1];        // first command line argument after name of program is put into element 1 of array: we use this for name of output file
        exponent = atoi(argv[2]);  //atoi: convert ascii input to integer
    }

    // Loop over powers of 10 to allow to rerun the calculations for different n values, n = # of mesh (aka integration) pts
    //So this way can right away perform same calculation with different # of mesh pts i.e. if i = 1, n = 10.  if i = 3, n=1000 etc.
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);         //pow(10.0) does 10^i and pow returns float by default. use pow(10.0,i) to return an int.
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i- when run the program (i.e. if typed 'test.exe out 10', it will make filenames be out1, out2 etc.
      fileout.append(argument);          //appends the 'argument' to our output file filename. The arguments will be i (i.e. 1, 2, 3)

      //   Define variables and pointers
      double h = 1.0/((double) n);  //define h=step size. 1.0 to make sure it returns a double. write: ((double) n) to make sure n is double also.
      double hh = h*h;              //define h^2 as hh
      // need to have dynamic memory allocation. For n points we have n+1 elements (arrays start from 0)
      double *d = new double[n+1];    //array  for diagonal matrix elements
      //double *e = new double[n+1];  //array for non-diagonal matrix elements, not needed for the special case
      double *f = new double[n+1];   //array for right hand side (then known fnc f(x))
      double *x = new double[n+1];   //array for x-values
      double *u = new double[n+1];   //array for solution

      //start clock timer
      high_resolution_clock::time_point start = high_resolution_clock::now();

      //Dirichlet boundary  conditions
      u[0] = u[n] = 0.0;
      //Special case, tridiagonal matrix
      d[0] = d[n] = 2.0;

      //Define x and f(x) values
      for (int i = 0; i <=n; i++){
          x[i] = i*h;
          f[i] = hh*ff(i*h);
      }
      //Precalculate
      for (int i = 1; i<n; i++) d[i] = (i+1.0)/((double) i);  //this is array d = (i+1)/i which saves FLOPs

      //Forward substition
      for (int i = 2; i <n; i++) f[i] = f[i] + f[i-1]/d[i-1];
      //for (int i = 2; i <n; i++) f[i] += f[i-1]/d[i-1];    // += is short way to write the above line and save floating pt operations

      //Backward substitution
      u[n-1] = f[n-1]/d[n-1]; //1st linear equation
      for (int i = n-2; i >0; i--) u[i] = (f[i]+u[i+1])/d[i];

      //stop clock timer and output time duration
      high_resolution_clock::time_point finish = high_resolution_clock::now();
      duration<double> time = duration_cast<duration<double>>(finish-start);
      cout << "time = " << time.count() << endl;

      //Setup output file
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase); //sets to write i.e. 10^6 as E6


      // Write to file:
      //     ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
         double RelativeError = fabs((exact(xval)-u[i])/exact(xval)); //fabs is absolute value
         ofile << setw(15) << setprecision(8) << xval;         //setprecision(8) sets that use 8 sigfigs
         ofile << setw(15) << setprecision(8) << u[i];
         ofile << setw(15) << setprecision(8) << exact(xval); //is analytical solution defined at top
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] d; delete [] f; delete [] u; //memory that was occupied by these arrays is now freed

    }

    return 0;   // returns zero upon success, this is not required but can be useful.
}


