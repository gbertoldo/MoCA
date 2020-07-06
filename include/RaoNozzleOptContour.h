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

#ifndef RAONOZZLEOPTCONTOUR_H
#define RAONOZZLEOPTCONTOUR_H

#include <vector>
#include "MoCToolBox.h"
// TODO (guilherme#1#): Verify RaoNozzleOptContour

/**

    Given the C- BD and the C+ DE (see 1958 Rao's paper), RaoNozzleOptContour

    provides the optimized contour BE of the nozzle.

*/
class RaoNozzleOptContour
{
    public:
        /**
            \brief Constructor
        */
        RaoNozzleOptContour(MoCToolBox& MoC, size_t Nit, double tol);


        /**
            \brief Destructor
        */
        virtual ~RaoNozzleOptContour();


        /**
            \brief Returns the MoCPoints on the Rao's optimized contour
        */
        std::vector<MoCPoint> optContour(const std::vector<MoCPoint>& CminusBD,
                                         const std::vector<MoCPoint>& CplusDE);

    private:

        /**
            \brief contourPoint calculates the next contour point, given a point p1 belonging to
            the contour and the points p2 and p3 that belong to C- or belong to a C+ (p2 and p3
            must be in the same characteristic). If the calculated point is between p2 and p3,
            isInside returns true, otherwise it returns false.
        */
        MoCPoint contourPoint(const MoCPoint& p1, const MoCPoint& p2, const MoCPoint& p3, bool& isInside);

    private:

        MoCToolBox            MoC; // MoC tool box
        size_t        Nit  {1000}; // Maximum number of iterations for iterative methods
        double        tol {1E-14}; // Tolerance for iteration errors

};

#endif // RAONOZZLEOPTCONTOUR_H
