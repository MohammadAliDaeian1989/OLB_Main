/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the computation of error norms.
 */


#include "olb3D.h"
#include "olb3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

using DESCRIPTOR = D3Q19<tag::MRT,FORCE>;


typedef enum {forced, nonForced} FlowType;

typedef enum {bounceBack, local, interpolated, bouzidi, freeSlip, partialSlip} BoundaryType;


// Parameters for the simulation setup
//FlowType flowType = nonForced;
BoundaryType boundaryType = interpolated;
int N = 25;            // resolution of the model
T Re = 200.;      // Reynolds number
T tau = 0.64265;
T tau_carr = 0.51;
T lx = 3.;             // channel lenght
T ly = 1.;             // channel width
T maxU = 0.0001;     // Max velocity
T Tmax = 2000;      // max. phys. time in smak
T Tprint = 1;      // Phys time at which the status of the system is print
// set the changes for n and m in powerLawBGKdynamics.h
T n = 0.15;              // parameter in power law model (n=1 Newtonian fluid)/ Or in the Carreau_Model
T n_carr=0.15;
T maxPhysT =1;
bool bcTypePeriodic = true; //true works only with one core

const T residuum = 1e-6;      // residuum for the convergence check

T mu_zero=1.0;
T mu_infinity=0;


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;
  T radius=1.0;
  T length=4.0;
  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, pipe);

    Vector<T, 3> origin(0, radius, radius);
    Vector<T, 3> extend = origin;

    // Set material number for inflow
    origin[0] = -converter.getPhysDeltaX() * 2;
    extend[0] = converter.getPhysDeltaX() * 2;
    IndicatorCylinder3D<T> inflow(origin, extend, radius);
    superGeometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = length - 2 * converter.getPhysDeltaX();
    extend[0] = length + 2 * converter.getPhysDeltaX();
    IndicatorCylinder3D<T> outflow(extend, origin, radius);
    superGeometry.rename(2, 4, 1, outflow);
  

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     CarreauUnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     Dynamics<T, DESCRIPTOR>& inDynamics,
                     Dynamics<T, DESCRIPTOR>& outDynamics,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(0), &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(1), &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(2), &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)

    sLattice.defineDynamics( superGeometry.getMaterialIndicator(3), &inDynamics );
 
    sLattice.defineDynamics( superGeometry.getMaterialIndicator(4), &outDynamics );
  
  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        CarreauUnitConverter<T,DESCRIPTOR> const& converter,
                        int iT, SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // Set initial and steady boundary conditions
  if ( iT==0 ) {

    // Define the analytical solutions for pressure and velocity
    T maxVelocity = converter.getCharLatticeVelocity();
    T distance2Wall = converter.getConversionFactorLength()/2.;
    //clout << "test11" <<std::endl;
    AnalyticalConst3D<T,T> rho( converter.getLatticeDensity( 1060 ));
      //  clout << "test22" <<std::endl;

    AnalyticalConst3D<T,T> u(0);
    //PowerLaw2D<T> u1( superGeometry, 3, 0, distance2Wall, ( n + 1. )/n );
   // clout << "test33" <<std::endl;

    // Set the analytical solutions for pressure and velocity
    AnalyticalConst3D<T,T> omega0( converter.getLatticeRelaxationFrequency() );
    //clout << "test1001" <<std::endl;
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 1, omega0 );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 3, omega0 );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 4, omega0 );
    //clout << "test44" <<std::endl;

    // Set the analytical solutions for pressure and velocity
    // Initialize all values of distribution functions to their local equilibrium

    sLattice.defineRhoU( superGeometry, 1, rho, u );
    sLattice.iniEquilibrium( superGeometry, 1, rho, u );

    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
    sLattice.defineRhoU( superGeometry, 3, rho, u );

    sLattice.iniEquilibrium( superGeometry, 4, rho, u );
    sLattice.defineRhoU( superGeometry, 4, rho, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}




// Output to console and files
void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer)
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "poiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );

    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );


  CarreauUnitConverter<T, DESCRIPTOR> const converter(
    int {N},            // resolution: number of voxels per charPhysL
    (T)   tau_carr,          // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   ly,           // charPhysLength: reference length of simulation geometry
    (T)   Re,           // Reynolds number
    (T)   n_carr,            // Carreauindex
    (T)   mu_zero,       // zero_viscosity
    (T)   1055.0           // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("poiseuille3d");


  // === 2nd Step: Prepare Geometry ===

  T radius=1;
  T length=4;
  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 2*singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 6;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getPhysDeltaX(), noOfCuboids);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );


// #if defined(ENABLE_MRT)
//   if (flowType == forced) {
//     bulkDynamics.reset(new ForcedMRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
//   }
//   else {
//     bulkDynamics.reset(new MRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
//   }
// #else
//   if (flowType == forced) {
//     bulkDynamics.reset(new ForcedBGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
//   }
//   else {
//     bulkDynamics.reset(new BGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
//   }
// #endif


  MRT_NN_dynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>(), mu_infinity, mu_zero, n_carr, converter.getConversionFactorLength(), converter.getConversionFactorTime());

  PeriodicPressureDynamics<T, DESCRIPTOR, MRT_NN_dynamics<T,DESCRIPTOR>> outDynamics( bulkDynamics,0.000000327,1,0);
  PeriodicPressureDynamics<T, DESCRIPTOR, MRT_NN_dynamics<T,DESCRIPTOR>> inDynamics( bulkDynamics,-0.000000327,-1,0);

  //prepareLattice and setBoundaryConditions
  prepareLattice(sLattice, converter, bulkDynamics,inDynamics,outDynamics, superGeometry);

  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  Vector<T, 3> center0Extended(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1Extended(length, radius, radius);


  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getConversionFactorLength()*N/11.);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss(sLattice, superGeometry, 2, converter, indicatorExtended);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, bulkDynamics, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
