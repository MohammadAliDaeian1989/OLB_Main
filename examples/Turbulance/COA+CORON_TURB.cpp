#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
//typedef D3Q19<> DESCRIPTOR;
using DESCRIPTOR = D3Q19<tag::MRT,FORCE>;



// Parameters for the simulation setup
const int N = 12;        // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 1600.; // max. simulation time in s, SI unit
const T radius=0.001; 


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );
  superGeometry.clean();
   Vector<T,3> extend;
  // extend={1.5,9.9,25.9};
  // Vector<T,3> origin;
  // origin={0.0063,6.4,25.0};
  // IndicatorCuboid3D<T> Outflow3(extend,origin);
  // superGeometry.rename(2,4,1,Outflow3);
  // IndicatorCircle3D<T> outflow4( 4.3,9.7,25.37,0,0,-1,0.7 );
  // IndicatorCylinder3D<T> layerOUTflow4( outflow4, 0.7);
  // superGeometry.rename( 2,4,1,layerOUTflow4);
  // IndicatorCircle3D<T> Inflow1(4.389,6.00019,11.2,0,0,1,1.60 );
  // IndicatorCylinder3D<T> layerinflow1( Inflow1, 0.7 );
  // superGeometry.rename( 2,3,1,layerinflow1);
  // IndicatorCircle3D<T> outflow2(4.760,10.0,0.0,0,0,1,1.0 );
  // IndicatorCylinder3D<T> layerOutflow2( outflow2, 0.9);
  // superGeometry.rename( 2,4,1,layerOutflow2);
  // IndicatorCircle3D<T> outflow5( 4.4,7.9,25.37, 0.,0.,-1., 0.33 );
  // IndicatorCylinder3D<T> layerOutflow5( outflow5, 2.0 );
  // superGeometry.rename( 2,4,1,layerOutflow5);
  // IndicatorCircle3D<T> outflow7(6.72026,7.62425,13.65,-1,0,0,0.4 );
  // IndicatorCylinder3D<T> layerOutflow7( outflow7,0.4);
  // superGeometry.rename( 2,5,1,layerOutflow7);
  // IndicatorCircle3D<T> outflow6(6.7,8.1,13.70,-1,0,0,0.2 );
  // IndicatorCylinder3D<T> layerOutflow6( outflow6,0.2 );
  // superGeometry.rename( 2,5,1,layerOutflow6);
  // clout << "material6" << std::endl;
  // IndicatorCircle3D<T> outflow8( 10.2345,3.537,11.6,0,1, 0,0.35 );
  // IndicatorCylinder3D<T> layerOutflow8( outflow8, 0.5 );
  // superGeometry.rename( 2,5,1,layerOutflow8 );
  // clout << "material7" << std::endl;
  // IndicatorCircle3D<T> outflow9( 9.6037,3.537,1.2274,0,1, 0,0.35 );
  // IndicatorCylinder3D<T> layerOutflow9( outflow9, 0.5 );
  // superGeometry.rename( 2,5,1,layerOutflow9 );
  // clout << "material8" << std::endl;
  // IndicatorCircle3D<T> outflow10( 4.46,8.977,8.291,-1,0,0, 0.3 );
  // IndicatorCylinder3D<T> layerOutflow10( outflow10, 0.5);
  // superGeometry.rename( 2,5,1,layerOutflow10 );
  // clout << "material9" << std::endl;
  // IndicatorCircle3D<T> outflow11( 4.111,5.79,6.87,0,1,0, 0.2 );
  // IndicatorCylinder3D<T> layerOutflow11( outflow11, 0.5 );
  // superGeometry.rename( 2,5,1,layerOutflow11 );
  // clout << "material10" << std::endl;

  extend={0.015,0.099,0.259};
  Vector<T,3> origin;
  origin={0.000063,0.064,0.250};
  IndicatorCuboid3D<T> Outflow3(extend,origin);
  superGeometry.rename(2,4,1,Outflow3);
  IndicatorCircle3D<T> outflow4( 0.043,0.097,0.2537,0,0,-1,0.007 );
  IndicatorCylinder3D<T> layerOUTflow4( outflow4, 0.007);
  superGeometry.rename( 2,4,1,layerOUTflow4);
  IndicatorCircle3D<T> Inflow1(0.04389,0.0600019,0.112,0,0,1,0.016 );
  IndicatorCylinder3D<T> layerinflow1( Inflow1, 0.007 );
  superGeometry.rename( 2,3,1,layerinflow1);
  IndicatorCircle3D<T> outflow2(0.04760,0.1,0.0,0,0,1,0.01 );
  IndicatorCylinder3D<T> layerOutflow2( outflow2, 0.009);
  superGeometry.rename( 2,4,1,layerOutflow2);
  IndicatorCircle3D<T> outflow5( 0.044,0.079,0.2537, 0.,0.,-1., 0.0033 );
  IndicatorCylinder3D<T> layerOutflow5( outflow5, 0.02 );
  superGeometry.rename( 2,4,1,layerOutflow5);
  IndicatorCircle3D<T> outflow7(0.0672026,0.0762425,0.1365,-1,0,0,0.004 );
  IndicatorCylinder3D<T> layerOutflow7( outflow7,0.004);
  superGeometry.rename( 2,5,1,layerOutflow7);
  IndicatorCircle3D<T> outflow6(0.067,0.081,0.1370,-1,0,0,0.002 );
  IndicatorCylinder3D<T> layerOutflow6( outflow6,0.002 );
  superGeometry.rename( 2,5,1,layerOutflow6);
  clout << "material6" << std::endl;
  IndicatorCircle3D<T> outflow8( 0.102345,0.03537,0.116,0,1, 0,0.0035 );
  IndicatorCylinder3D<T> layerOutflow8( outflow8, 0.005 );
  superGeometry.rename( 2,5,1,layerOutflow8 );
  clout << "material7" << std::endl;
  IndicatorCircle3D<T> outflow9( 0.096037,0.03537,0.12274,0,1, 0,0.0035 );
  IndicatorCylinder3D<T> layerOutflow9( outflow9, 0.005 );
  superGeometry.rename( 2,5,1,layerOutflow9 );
  clout << "material8" << std::endl;
  IndicatorCircle3D<T> outflow10( 0.0446,0.08977,0.08291,-1,0,0, 0.003 );
  IndicatorCylinder3D<T> layerOutflow10( outflow10, 0.005);
  superGeometry.rename( 2,5,1,layerOutflow10 );
  clout << "material9" << std::endl;
  IndicatorCircle3D<T> outflow11( 0.04111,0.0579,0.0687,0,1,0, 0.002 );
  IndicatorCylinder3D<T> layerOutflow11( outflow11, 0.005 );
  superGeometry.rename( 2,5,1,layerOutflow11 );
  clout << "material10" << std::endl;



  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
  //////auto bulkIndicator2 = superGeometry.getMaterialIndicator({3});
  //////auto bulkIndicator3 = superGeometry.getMaterialIndicator({4});
  ////// auto bulkIndicator4 = superGeometry.getMaterialIndicator({2});

  sLattice.defineDynamics( bulkIndicator, &bulkDynamics );
  ////sLattice.defineDynamics( bulkIndicator2, &bulkDynamics );
  ////sLattice.defineDynamics( bulkIndicator3, &bulkDynamics );
  ////sLattice.defineDynamics( bulkIndicator4, &bulkDynamics );


  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

    setBouzidiZeroVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2, stlReader);

  // Setting of the boundary conditions

  //// if local boundary conditions are chosen
 	////setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
 	////setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

 	//if interpolated boundary conditions are chosen
 	setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
 	setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);
	setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF(1);
 
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );

  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  ////int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTmaxStart = 40;
  int iTupdate = 5;
  clout << "BoundaryValue ... check1" << std::endl;
  //if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] =converter.getCharLatticeVelocity();


  AnalyticalConst3D<T,T> U_Bound( T( 0 ), T( 0 ), maxVelocity[0] );
  AnalyticalConst3D<T,T> P_Bound( 1.0 );

  sLattice.defineU  ( superGeometry, 3,U_Bound  );
  sLattice.defineRho( superGeometry, 4,P_Bound  );
  sLattice.defineRho( superGeometry, 5,P_Bound  );
  

  clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
  //}
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer,
                 STLreader<T>& stlReader )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "yshapeduct" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 5 );
  SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR> quality( sLattice, converter );
  SuperRoundingF3D<T, T> roundedQuality ( quality, RoundingMode::NearestInteger );
  SuperDiscretizationF3D<T> discretization ( roundedQuality, 0., 2. );

  clout << "getresult ... check1" << std::endl;

  vtmWriter.addFunctor( quality );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( yPlus );
  clout << "getresult ... check2" << std::endl;
  //const int vtkIter  = converter.getLatticeTime( .3 );
  const int vtkIter  =500;
  //const int statIter = converter.getLatticeTime( .1 );
  clout << "getresult ... check3" << std::endl;
  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
clout << "getresult ... check4" << std::endl;
    vtmWriter.createMasterFile();
  }
clout << "getresult ... check5" << std::endl;
  // Writes the vtk files
  if ( iT%vtkIter == 0 ) {
    vtmWriter.write( iT );

    {
      SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
       //write output as JPEG
      heatmap::write(planeReduction, iT);
    }
clout << "getresult ... check5" << std::endl;
    {
      BlockReduction3D2D<T> planeReduction( discretization, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> jpeg_scale;
      jpeg_scale.colour = "blackbody";
      jpeg_scale.name = "quality";
      heatmap::write( planeReduction, iT, jpeg_scale );
    }
  }
clout << "getresult ... check5" << std::endl;
//////  // Writes output on the console
//////  ////if ( iT%statIter == 0 ) {
//////    // Timer console output
//////    timer.update( iT );
//////    timer.printStep();
//////clout << "getresult ... check6" << std::endl;
//////    // Lattice statistics console output
//////    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
//////
//////    // Drag, lift, pressure drop
//////    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );
//////    SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );
//////clout << "getresult ... check7" << std::endl;
//////    std::vector<T> point1V = superGeometry.getStatistics().getCenterPhysR( 5 );
//////    std::vector<T> point2V = superGeometry.getStatistics().getCenterPhysR( 5 );
//////    T point1[3] = {};
//////    T point2[3] = {};
//////    for ( int i = 0; i<3; i++ ) {
//////      point1[i] = point1V[i];
//////      point2[i] = point2V[i];
//////    }
//////    point1[0] = superGeometry.getStatistics().getMinPhysR( 5 )[0] - converter.getConversionFactorLength();
//////    point2[0] = superGeometry.getStatistics().getMaxPhysR( 5 )[0] + converter.getConversionFactorLength();
//////clout << "getresult ... check8" << std::endl;
//////    T p1, p2;
//////    intpolatePressure( &p1,point1 );
//////    intpolatePressure( &p2,point2 );
//////
//////    clout << "pressure1=" << p1;
//////    clout << "; pressure2=" << p2;
//////
//////    T pressureDrop = p1-p2;
//////    clout << "; pressureDrop=" << pressureDrop;
//////
//////    T dragA[3];
//////    int input1[0];
//////    drag( dragA, input1 );
//////    clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << endl;
//////
//////    int input[4] = {};
//////    SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
//////    T yPlusMax[1];
//////    yPlusMaxF( yPlusMax,input );
//////    clout << "yPlusMax=" << yPlusMax[0] << endl;
//////clout << "getresult ... check9" << std::endl;
//////  }
}

int main( int argc, char* argv[] )
{
  ofstream myfile;
  myfile.open("Time.txt") ;
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

 ///* if (argc > 1) {
 //   if (argv[1][0]=='-'&&argv[1][1]=='h') {
 //     OstreamManager clout( std::cout,"help" );
 //     clout<<"Usage: program [Resolution] [FlowType] [BoundaryType]"<<std::endl;
 //     clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
 //     clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated, 3=bouzidi, 4=freeSlip, 5=partialSlip"<<std::endl;
 //     clout<<"Default: FlowType=forced, Resolution=21, BoundaryType=bouzidi"<<std::endl;
 //     return 0;
 //   }
 // }

 // if (argc > 1) {
 //   N = atoi(argv[1]);
 //   if (N < 1) {
 //     std::cerr << "Fluid domain is too small" << std::endl;
 //     return 1;
 //   }
 // }

 // if (argc > 2) {
 //   int flowTypeNumber = atoi(argv[2]);
 //   if (flowTypeNumber < 0 || flowTypeNumber > (int)nonForced) {
 //     std::cerr << "Unknown fluid flow type" << std::endl;
 //     return 2;
 //   }
 //   flowType = (FlowType) flowTypeNumber;
 // }

 // if (argc > 3) {
 //   int boundaryTypeNumber = atoi(argv[3]);
 //   if (boundaryTypeNumber < 0 || boundaryTypeNumber > (int) partialSlip) {
 //     std::cerr << "Unknown boundary type" << std::endl;
 //     return 3;
 //   }
 //   boundaryType = (BoundaryType) boundaryTypeNumber;
 // }*/

//  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
//     int {N},                // resolution: number of voxels per charPhysL
//     (T)   0.51,             // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
//     (T)   0.001,            // charPhysLength: reference length of simulation geometry
//     (T)   0.1,               // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
//     (T)   0.003/1050,   // physViscosity: physical kinematic viscosity in __m^2 / s__
//     (T)   1050.0               // physDensity: physical density in __kg / m^3__
//   );

  UnitConverter<T,DESCRIPTOR> converter(
    (T)   0.001/N,        // physDeltaX: spacing between two lattice cells in __m__
    (T)   0.000012074/9,     // physDeltaT: time step in __s__
    (T)   0.001,       // charPhysLength: reference length of simulation geometry 
    (T)   0.6,         // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.003/1055.,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1055           // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  clout << "Conversionlength  " <<  converter.getConversionFactorLength()    << std::endl;
 converter.print();
  // Writes the converter log in a file
  converter.write("COA-COR_71");
  clout << "Check11 ..." << std::endl;
  STLreader<T> stlReader( "COA-COR_71.stl", converter.getConversionFactorLength(),0.001);
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  clout << "Check22 ..." << std::endl;

   const int noOfCuboids = 16;

  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );
   // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

    //BGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );

    //MRTdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(),
      instances::getBulkMomenta<T, DESCRIPTOR>(), 0.1 );

  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
   timer.start();
  //prepareLattice and set boundaryCondition
   prepareLattice( sLattice, converter, bulkDynamics, stlReader, superGeometry );


    for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {


  auto start = std::chrono::system_clock::now();
   setBoundaryValues( sLattice, converter, iT, superGeometry );
   sLattice.collideAndStream();
   auto end = std::chrono::system_clock::now();
   getResults( sLattice, converter,iT, superGeometry, timer, stlReader );
	auto end_write = std::chrono::system_clock::now();
	auto elapsed = end - start;
	auto elapsed1 = end_write - start;
    //std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
      //        << 1000.0 * (c_end - c_start);
   myfile << iT << ","  <<  elapsed.count() << ","  << elapsed1.count()<< endl;

	}

  return 0;
 myfile.close();
 timer.stop();


}
