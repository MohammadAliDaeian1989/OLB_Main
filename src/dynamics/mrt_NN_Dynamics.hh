/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef MRT_DYNAMICS_NN_HH
#define MRT_DYNAMICS_NN_HH

#include <algorithm>
#include <limits>
#include "mrtHelpers.h"

namespace olb {

//==============================================================================//
/////////////////////////// Class MRT_NN_dynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param lambda_ will be used as an
 */

// Original implementation based on:
// D'Humieres et al., "Multiple-relaxation-time lattice Boltzmann models in three dimensions",
// Phil: Trans. R. soc. Lond. A (2002) 360, 437-451
// and
// Yu et al,, "LES of turbulent square jet flow using an MRT lattice Boltzmann model",
// Computers & Fluids 35 (2006), 957-965
template<typename T, typename DESCRIPTOR>
MRT_NN_dynamics<T,DESCRIPTOR>::MRT_NN_dynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ ,T mu_infinity1, T mu_zero1, T n_carr1, T dx1, T dt1)
  : BasicDynamics<T,DESCRIPTOR>(momenta_), omega(omega_), lambda(omega_)
{

    mu_infinity=mu_infinity1;
    mu_zero=mu_zero1;
    n_carr=n_carr1;
    dx=dx1;
    dt=dt1;
  this->getName() = "MRT_NN_dynamics";  
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s<T,DESCRIPTOR>(iPop);
  }
  for (int iPop  = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = 1;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) *
                                rt[kPop];
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
T MRT_NN_dynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void MRT_NN_dynamics<T,DESCRIPTOR>::computeAllEquilibrium(T momentaEq[DESCRIPTOR::q],
                                                T rho, const T u[DESCRIPTOR::d],
                                                const T uSqr)
{
  mrtHelpers<T,DESCRIPTOR>::computeEquilibrium(momentaEq, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void MRT_NN_dynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{ OstreamManager clout(std::cout, "MRT_NN");
  typedef DESCRIPTOR L;
  typedef mrtHelpers<T,DESCRIPTOR> mrtH;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  const auto oldOmega = cell.template getField<descriptors::OMEGA>();
  const auto newOmega = this->computeOmegaCarr(cell, oldOmega, rho, mu_infinity, mu_zero, n_carr, dx, dt);
  cell.template setField<descriptors::OMEGA>(newOmega);
  //clout << "OMEGA_NEW" << "," <<newOmega <<"," << 1/newOmega << std::endl;
  for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
    for (int iPop  = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
      //clout << "========>" << descriptors::shearIndexes<DESCRIPTOR>() << std::endl;
      //clout << descriptors::shearViscIndexes<DESCRIPTOR>(iPop)<< std::endl;
      //clout << "1>>"<<invM_S[kPop][descriptors::shearViscIndexes<DESCRIPTOR>(iPop)]<< std::endl;
      invM_S[kPop][descriptors::shearViscIndexes<DESCRIPTOR>(iPop)]=descriptors::invM<T,DESCRIPTOR>(kPop,descriptors::shearViscIndexes<DESCRIPTOR>(iPop))*newOmega;
      //clout << "2>>"<< invM_S[kPop][descriptors::shearViscIndexes<DESCRIPTOR>(iPop)]<< std::endl;
   }
  }
  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T MRT_NN_dynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void MRT_NN_dynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR>
T MRT_NN_dynamics<T,DESCRIPTOR>::getLambda() const
{
  return lambda;
}

template<typename T, typename DESCRIPTOR>
void MRT_NN_dynamics<T,DESCRIPTOR>::setLambda(T lambda_)
{
  lambda = lambda_;
}



template<typename T, typename DESCRIPTOR>
T MRT_NN_dynamics<T,DESCRIPTOR>::computeOmegaCarr(Cell<T,DESCRIPTOR>& cell, T omega0,
           T rho,T mu_ifinity,T mu_zero,T n_carr, T dx, T dt)           
{

  T a_carr=0.644;
  T Lanndha_carr=0.11*dt;

  T pre2 = pow(descriptors::invCs2<T,DESCRIPTOR>()/2.* omega0/rho,2.); // strain rate tensor prefactor
  //T gamma = sqrt(2.*pre2*PiNeqNormSqr(cell)); // shear rate
  T gamma = sqrt(2.*pre2*lbHelpers<T,DESCRIPTOR>::computePiNeqNormSqr(cell)); // shear rate
  //std::cout<< "-111111111->"  <<  gamma << std::endl;
  //std::cout<< "-------->"  <<  Lanndha_carr << std::endl;
  // T nuNew = _m*pow(gamma,_n-1.); //nu for non-Newtonian fluid
  //T nuNew = (numin_carr-numix_carr)*(1+pow((Landha*gamma)^(a_carr))^((n_carr-1)/(a_carr)); //nu for non-Newtonian fluid

  T nuNew1 = ((mu_zero-mu_ifinity)*pow((1+pow(gamma*Lanndha_carr,a_carr)),((n_carr-1)/(a_carr)))+mu_ifinity);
  T nuNew =nuNew1/1055.0* dt/dx/dx;
    
  T newOmega = 1./(nuNew*descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
  //std::cout<< "-------->"  <<  dt << "," << dx << std::endl;
  //std::cout<< "-------->"  << nuNew1<<"-------->"<< nuNew << "--------> " << newOmega << "-----> " << Lanndha_carr<< "      ";
  
  return newOmega;
  //return omega0;
}




template<typename T, typename DESCRIPTOR>
MRT_NN_dynamics2<T,DESCRIPTOR>::MRT_NN_dynamics2 (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : MRT_NN_dynamics<T,DESCRIPTOR>(omega_, momenta_)
{
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s_2<T,DESCRIPTOR>(iPop);
  }
  for (int iPop  = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S_2[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S_2[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) *
                                  rt[kPop];
        }
      }
    }
  }
}

// Stabalized MRT scheme with uniform relaxation times
template<typename T, typename DESCRIPTOR>
void MRT_NN_dynamics2<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef mrtHelpers<T,DESCRIPTOR> mrtH;

  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S_2);

  statistics.incrementStats(rho, uSqr);
}



template<typename T, typename DESCRIPTOR>
ForcedMRT_NN_dynamics<T,DESCRIPTOR>::ForcedMRT_NN_dynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : MRT_NN_dynamics<T,DESCRIPTOR>(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR>
void ForcedMRT_NN_dynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtHelpers<T,DESCRIPTOR>::mrtCollision(cell, rho, u, this->invM_S);
  mrtHelpers<T,DESCRIPTOR>::addExternalForce(cell, rho, u, this->invM_S);

  statistics.incrementStats(rho, uSqr);
}


} // end namespace

#endif

