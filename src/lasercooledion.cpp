/**
 * @file lasercooledion.cpp
 * @brief Function definitions for class derived from Ion
 *
 * @author Chris Rennick
 * @copyright Copyright 2014 University of Oxford.
 */

/**
 *  @class LaserCooledIon
 *  @brief Representation of a trapped and laser cooled ion.
 *  This class extends the trapped ion class to include frictional cooling
 *  forces and photon recoil heating force.
 *
 *  Heating term arising from photon recoil is implemented as a Langevin process
 *  with a Gaussian momentum distribution.
 */

#include "include/ion.h"

#include "include/ccmdsim.h"

#include <math.h>

#include <stdlib.h>

#include <cmath>

#include "include/vector3D.h"

using namespace std;

/**
 *  Construct a new laser cooled ion. The trap  are passed up to the
 *  `TrappedIon` parent class, and laser cooling parameters are stored.
 *  @param ion_trap A pointer to the ion trap.
 *  @param type     A pointer to ion parameters.
 */
LaserCooledIon::LaserCooledIon(const IonTrap_ptr ion_trap, const IonType& type,
        const SimParams& sp)
    : TrappedIon(ion_trap, type), heater_(sp.random_seed) {
    heater_.set_kick_size(sqrt(ionType_.recoil));
}

/**
 *  @brief Change the ion velocity due to the laser cooling and trapping forces.
 *  The trapping force is handled first by calling the parent class kick
 *  function. The radiation pressure and friction force (cooling) are calculated
 *  and applied.
 *
 *  @param dt   Time step.
 */
inline void LaserCooledIon::kick(double dt) {
    this->TrappedIon::kick(dt);

    // 1D radiation pressure force.
    Vector3D pressure(0.0, 0.0, 0.015);
    // Randomly apply pressure from positive or negative \c z direction,
    // weighted by the `direction` parameter.
    if (heater_.kick_direction(ionType_.direction))
        this->Ion::kick(dt, pressure);
    else
        this->Ion::kick(dt, -pressure);

    // 1D Laser cooling friction force
    // This force must be evaluated last to allow its effect to be
    // undone by the call to velocity_scale
    Vector3D friction = get_friction();
    this->Ion::kick(dt, -friction);

    return;
}

//added function
float fscatt(Vector3D vel, const LaserParams& lp, const IonType& type) {
    
    float pi = 3.14159265359;
    double Gamma = type.A21;
    double IdIsat = lp.IdIsat;
    double delta = lp.delta;
	double Gamma3 = pow(Gamma,3);
	double Gamma2 = std::pow(Gamma,2);
	Vector3D k = (0,0,(2*pi) / lp.wavelength );
    
    double gamma = 0.5 * (Gamma3);
    gamma *= IdIsat;
    double x = delta - dot(vel,k);
    double x2 = std::pow(x,2);
    gamma /= (Gamma2 + (4 * x2));
    return gamma;
	
}

// added function
void isoEmit(Vector3D vel, const LaserParams& lp){
    
    float h = 6.62607*std::pow(10,-34);
    float pi = 3.14159265359;
    
    float rnumu = random.random();
    float rnumv = random.random();
    float theta = 2 * pi * rnumu;
    float phi = acos((2 * rnumv) -1);
    float x = cos(theta) * sin(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(phi);
    
    vel -= (x * (h/lp.wavelength), y * (h/lp.wavelength), z * (h/lp.wavelength));
    
    return;

}

//added function
void Emit(Vector3D vel, const IonType& type, const LaserParams& lp) {
    
    float dt = 1.0*std::pow(10,-9);
    float fs = fscatt(vel, lp, type); //Probability of stimulated emission s^-1
    fs += 1.4*std::pow(10,-8); //Probability of spontaneous emission s^-1 from NIST
    fs *= dt;
    
    //number = sum(1 for item in r if item <= fscatt);
    
    //for i in range(number){
    isoEmit(vel, lp);
    return;
}	

// added function  
void Absorb(Vector3D vel, const IonType& type, const LaserParams& lp){

    float h = 6.62607*std::pow(10,-34);
    float dt = 1.0*std::pow(10,-9);
    float fs = fscatt(vel, lp, type); //Probability of stimulated absorption s^-1
    fs *= dt;
    
    //double number = sum(1 for item in r if item <= fscatt);
    
    float recoil_momentum = (h/lp.wavelength);// * number
    vel -= (recoil_momentum,0,0);

    return;
	
}

/**
 *  @brief Calculate laser cooling friction vector.
 *
 *  @return Friction vector.
 */
Vector3D LaserCooledIon::get_friction() const {
    return Vector3D(0.0, 0.0, ionType_.mass*ionType_.beta*vel_.z);
}

/**
 *  @brief Correct for friction forces in Velocity Verlet algorithm.
 *  This follows the algorithm given by:
 *      M. Tuckerman and B. J. Berne,
 *      J. Chem. Phys. 95, 4389 (1991), Eqn.
 *  Note that this routine is called with dt/2, so use of "dt" referring to the
 *  complete timestep requires multiplication by 2.
 *
 *  @param dt   Time step.
 */

void LaserCooledIon::velocity_scale(double dt) {
    // Note that this routine is called with dt/2, so use
    // of "dt" referring to the complete timestep requires
    // multiplication by 2.
    // double two_dt = 2.0*dt;

    // Eqn. 3.6, undo friction term
    Vector3D friction = get_friction();
    this->Ion::kick(dt, friction);

    // Eqn. 3.7
    // vel_.z /= 1.0 + two_dt*dt*ionType_.beta/ionType_.mass;
    vel_.z /= 1.0 + dt*ionType_.beta/ionType_.mass;
}

/**
 *  @brief Change velocity due to heating.
 *  The photon recoil force is generated as a random vector.
 *
 *  @param dt   Time step size.
 */
void LaserCooledIon::heat(double dt) {
    Vector3D heating = this->heater_.random_kick();
    this->Ion::kick(dt, heating);
}


