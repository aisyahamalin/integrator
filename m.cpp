//
//  main.cpp
//  frog
//
//  Created by Aisyah Zainal on 05/07/2020.
//  Copyright © 2020 Aisyah Zainal. All rights reserved.
//

//ofstream: outputfilestream: to create and write to files
//ifstream: inputfilestream: to read from files
//fstream: filestream: to both write to and read from files

#include<iostream>  //for cout, cin
#include<fstream>   //for writing and reding from files
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

//a function to make dimensionless

//a function to convert back to actual units





//=========================================================================
//specifying the potential class incorporating both getpot and getforce
//this will give the potential and force
class potential {
private:
    ;
public:
    double getpot(double *pos);//function inside class//*pos =>> declaring a pointer inside brackets//pos is position inside array
    void getforce(double *pos, double *res);   //
};
//--------------------------------------------------------------------------
//a function returning the Kepler potential, depending on r
double potential::getpot(double *pos) {
        double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        return (-GM/(r));
}
//--------------------------------------------------------------------------
//a function giving the force, for each dimension x, y, z
//remember that void function doesn't return a value
void potential::getforce(double *pos, double *force) {
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double rm3 = 1.0/(pow(r, 3.0));

    for (int i = 0; i < DIMENSION; i++) {
        force[i] = -GM*pos[i]*rm3;
    }
}
//======================================================================
//========================================================================
//specifying the star class incoporating the functions
//The object of mass M that effect the orbit of point mass m
class star {
private:
    double q[DIMENSION]; //an array with 3 elements (position)
    double p[DIMENSION]; //an array with 3 elements (velocity)
public:
    void setstar(double *x, double *v);
    void printcoords();
    
    void drift(double h);
    void kick(double h, double *force);
    void getforce(double *force, potential *Phi);
    void leapfrog(double h, potential *Phi);
    
    void printfile(ofstream& fileout);
    
    //defining the function for varying h
    double settime(double e);
    
    //function for evaluating the energy
    //double getE();
    
    //function for evaluating the angular momentum
    
    //function for applying the RK4 method
    //void runge-kutta();
    
};
//-----------------------------------------------------------------------
/*
//calculating the kinetic energy and potential, to get total energy---------------------------------------------------------
double star::getE(potential *Phi) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        Ekin += 0.5*p[i]*p[i]; //0.5*velocity^2
    }

    return (Ekin + Phi->getpot(q));
}//----------------------------------------------------------------------------------------------------
*/
    


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

//void function for the drift [updating the position of x, y, z each time]
void star::drift(double h) {
        for (int i = 0; i < DIMENSION; i++) {
            q[i] += (h*p[i]);
        } //position += h*velocity
}
//void function for the kick [updating the velocity using force/acceleration]
void star::kick(double h, double *force) {
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += h*force[i];
    } //velocity += h*force  (force is the acceleration with point mass m)
}
//void function for the leapfrog method!!!
void star::leapfrog(double h, potential *Phi) {
    double h2 = h*0.5;
    double *force = new double[3];  //'new' operator requesting to allocate memory dynamically, here there are 3 elements...
    getforce(force, Phi);
    drift(h2);        //drift
    kick(h, force);   //kick
    drift(h2);        //drift
    delete [] force; //deleting memory so it won't take up space allocated by the 'new' operator
}


void star::printfile(ofstream& fileout){
  fileout << q[0] << " " << q[1] << " " << q[2] << " " << p[0] << " " << p[1] << " " << p[2]  << endl;
}
//void function for
void star::getforce(double *force, potential *Phi) {
    Phi->getforce(q, force);  //calling getforce with position q and force
}                             // -> is the arrow operator, to reference individual members of classes
//here, we access a member of the class through a pointer
//Phi is a pointer
//a pointer stores the address of a variable
// writing *Phi is obtaining the value of variable, not address
//=========================================================================



//THE MAIN FUNCTION
int main() {
    double *x = new double[DIMENSION];  //value of pointer allocated dynamical memory
    double *v = new double[DIMENSION];  //value of pointer allocated dynamical memory
    //the values of each array
    x[0] = 8.0; x[1] = 0.0; x[2] = 0.0;
    v[0] = 0.0; v[1] = 1.0; v[2] = 0.3;

    star pedro; //calling it pedro
    pedro.setstar(x, v); //setting the position x and velocity v for pedro

    delete [] x; //deleting the allocated memory
    delete [] v; //deleting the allocated memory
    pedro.printcoords();

    //defining the step-length //which is a constant throughout code
    double h = 1.0e-6;  //unit of timestep is: 1 Myr
    //here, h is 1 year
    
    //defining the parameter constant e
    double e=10e-6;
    double dt;

    potential Phi; //calling it Phi //making an instance of potential called phi

    //opening the file to print to
    ofstream myfile ("coor.dat", std::ios_base::app);
    for (int i = 0; i < 1000; i++) {   //one thousand times
        for (int ii = 0; ii < 300000; ii++) { //for each one time out of a thousand, do it 300000 times
            dt = pedro.settime(e);
                
            pedro.leapfrog(dt, &Phi); //the leapfrog for that particular star
        }                   //address of variable's value: &Phi

        //printing the new coordinates after applying the leapfrog
        pedro.printcoords();
        pedro.printfile(myfile);
            }
    myfile.close();
    return 0;
}
