
#include<iostream>  //for cout, cin
#include<fstream>   //for writing and reading from files
#include<string>
#include<functional>
#include<cmath>
#include<random>
#include<thread>
#include<mutex>
#include<algorithm>

#include <math.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass

using namespace std;
using std::setprecision;  //for decimal precision...

//all about units
//so that the units are dimensionless in the programming
//M_u = 5.6*e10 Solar_mass
//R_u = 3.5 kpc
//V_u = 262 km/s
//T_u = 0.013 Gyr


//=========================================================================
//specifying the potential class incorporating both getpot and getforce
//this will give the potential and force
class potential {
private:
    ;
public:
    double getpot(double *pos); //a getter function
    void getforce(double *pos, double *res); //a getter function
};
//--------------------------------------------------------------------------
//defining the getter functions
//a function returning the Kepler potential, depending on r
double potential::getpot(double *pos) {
        double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        return (-GM/(r));
}
//--------------------------------------------------------------------------
//a function giving the force, for each dimension x, y, z
void potential::getforce(double *pos, double *force) {
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double rm3 = 1.0/(pow(r, 3.0));

    for (int i = 0; i < DIMENSION; i++) {
        force[i] = -GM*pos[i]*rm3;
    }
}
//======================================================================
//========================================================================
class star {
private:
    double q[DIMENSION]; //an array with 3 elements (position)
    double p[DIMENSION]; //an array with 3 elements (velocity)
    double E;
    double L[DIMENSION]; //an array for the angular momentum, in each direction
    double k1r[DIMENSION], k2r[DIMENSION],k3r[DIMENSION], k4r[DIMENSION];
    double k1v[DIMENSION], k2v[DIMENSION],k3v[DIMENSION], k4v[DIMENSION];

public:
    void setstar(double *x, double *v); //a constructor
    void printcoords(); //a constructor
    void getforce(double *force, potential *Phi);
    void printfile(ofstream& fileout);
    //defining the function for varying h
    double settime(double e);
    //function for evaluating the energy
    double getE(potential *Phi);
    //function for evaluating the angular momentum
    void getL(); //position and velocity are already private!
    
    //function for applying the RK4 method
    void runge_kutta(double h, potential *Phi);
    //RK_stepA
    void stepA(double h, double *force);
    void stepB(double h, double *force);
    void stepC(double h, double *force);
    void stepD(double h, double *force);
};
//-----------------------------------------------------------------------

//===========================================================
//RK4 method
void star::runge_kutta(double h, potential *Phi){
    
    double *force = new double[3]; //creating memory from the pointer
    getforce(force, Phi);
    stepA(h,force);
    stepB(h,force);
    stepC(h,force);
    stepD(h,force);
    
    for (int i = 0; i < DIMENSION; i++){
        //position q update
        q[i] += h/6.0 *(k1r[i] + k2r[i] + k3r[i] + k4r[i]);
        //velocity p update
        p[i] += h/6.0 *(k1v[i] + k2v[i] + k3v[i] + k4v[i]);
    }
    
    E = getE(Phi);
    //E*=1.0e11;
    //L = getL();
    
    delete [] force;
    
}
//====================================================================
//in each step, loop over each array element
//the derivative in the k1 direction
void star::stepA(double h, double *force){
    for (int i = 0; i < DIMENSION; i++){
        k1r[i] = p[i];     //k1r = velocity
        k1v[i] = force[i];     //k1v = acceleration
    }
}
//the derivative in the k2 direction
void star::stepB(double h, double *force){
    for (int i = 0; i < DIMENSION; i++){
        k2r[i] = p[i] + k1r[i]*h/2.0;   //velocity + k1r h/2
        k2v[i] = force[i] + k1v[i]*h/2.0;   //accel + k1v h/2
    }
}
//the derivative in the k3 direction
void star::stepC(double h, double *force){
    for (int i = 0; i < DIMENSION; i++){
        k3r[i] = p[i] + k2r[i]*h/2.0;//velocity + k2r h/2
        k3v[i] = force[i] + k2v[i]*h/2.0;  //accel + k2v h/2
    }
}
//the derivative in the k4 direction
void star::stepD(double h, double *force){
    for (int i = 0; i < DIMENSION; i++){
        k4r[i] = p[i] + k3r[i]*h;     //velocity + k3r h
        k4v[i] = force[i] + k3v[i]*h;     //accel + k3v h
    }
}
//====================================================================

//angular momentum L
void star::getL(){ //args are position and velocity
    //double L[DIMENSION] = [0.0, 0.0, 0.0]; //to initialise L
    for (int i = 0; i < DIMENSION; i++) {
        L[i] += p[i]*q[i]; //velocity*position
    }
}

//calculating the kinetic energy and potential, to get total energy---------------------------------------------
double star::getE(potential *Phi) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        Ekin += 0.5*p[i]*p[i]; //0.5*velocity^2
    }
    return (Ekin + Phi->getpot(q));
}//----------------------------------------------------------------------------------------------------
//function for varying the time-step h, depending on r and v values of particle
double star::settime(double e){
    double mod_R = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    double mod_V = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    double dt;
    return dt = e*mod_R/mod_V;
}
//void function for assigning the position x and velocity v
void star::setstar(double *x, double *v) {
    for (int i = 0; i < DIMENSION; i++) {
        q[i] = x[i];   //assigning the position
        p[i] = v[i];   // assigning the velocity
    }
}
//void function for printing the coordinates to screen
void star::printcoords() {
    for (int i = 0; i < DIMENSION; i++) {
        cout << q[i] << " " << p[i] << " ";
    }   //the first column gives the position, the second one gives velocity
    cout << "\n";
}
//void function for
void star::getforce(double *force, potential *Phi) {
    Phi->getforce(q, force);  //calling getforce with position q and force
}
void star::printfile(ofstream& fileout){
  fileout << q[0] << " " << q[1] << " " << q[2] << " " << p[0] << " " << p[1] << " " << p[2]  << " " << E << " " << L[0] << " " << L[1] << " " << L[2] << endl;
}
//=========================================================================


//THE MAIN FUNCTION
int main() {
    double *x = new double[DIMENSION];  //value of pointer allocated dynamical memory
    double *v = new double[DIMENSION];  //value of pointer allocated dynamical memory
    //the values of each element
    x[0] = 8.0; x[1] = 0.4; x[2] = 0.0;
    v[0] = 0.1; v[1] = 1.0; v[2] = 0.2;

    star pedro; //calling it pedro
    pedro.setstar(x, v); //setting the position x and velocity v for pedro

    delete [] x; //deleting the allocated memory
    delete [] v; //deleting the allocated memory
    pedro.printcoords();

    //defining the step-length //which is a constant throughout code
    double h = 1.0e-6;  //unit of timestep is: 1 Myr //here, h is 1 year
    //my original 1.0e-6
    
    //defining the parameter constant e
    double e = 1.0e-6; //my original 1.0e-7
    double dt;

    potential Phi; //calling it Phi //making an instance of potential called phi

    //opening the file to print to
    ofstream myfile ("coor_RK.dat", std::ios_base::app);
    for (int i = 0; i < 1000; i++) {   //one thousand times //how many times you print to file
        for (int ii = 0; ii < 10000; ii++) { //for each one time out of a thousand, do it 300000 times
                                             //actually implementing the integrator
            dt = pedro.settime(e);
            pedro.runge_kutta(dt, &Phi); //applying the integrator
        }                             //
        
        //printing the new coordinates after applying the integrator
        pedro.printcoords();         
        pedro.printfile(myfile);
            }
    myfile.close();

    return 0;
}








