/**
    MoCA - Method of Characteristics for Axisymmetric flows

    Copyright (C) 2020 by

    \author: Guilherme Bertoldo
             Federal University of Technology - Parana (UTFPR)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "IsentropicFlowRelations.h"

#include <cmath>

// Constructor
IsentropicFlowRelations::IsentropicFlowRelations(const double & gamma)
{
   this->GM = gamma;
}


// Destructor
IsentropicFlowRelations::~IsentropicFlowRelations()
{
   //dtor
}

// Pressure ratio p/p0 (0 index means stagnation condition)
double IsentropicFlowRelations::Pr(const double& M){

   return pow(1.0+(GM-1.0)/2.0*M*M, -GM/(GM-1.0));

}

// Temperature ratio T/T0 (0 index means stagnation condition)
double IsentropicFlowRelations::Tr(const double& M){

   return 1.0/(1.0+(GM-1.0)/2.0*M*M);

}

// Density ratio rho/rho0 (0 index means stagnation condition)
double IsentropicFlowRelations::rhor(const double& M){

   return pow(1.0+(GM-1.0)/2.0*M*M, -1.0/(GM-1.0));

}
