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

#include "RaoControlSurface.h"
#include <algorithm>
#include "NumRootFinding.h"
#include "NumIntegration.h"

RaoControlSurface::RaoControlSurface(const double& gamma,
                                     const double& Par,
                                     const double& ME,
                                     MoCToolBox& MoC,
                                     const size_t N):
    gamma(gamma),
    Par(Par),
    N(N),
    ME(ME),
    MoC(MoC),
    isenFlow(gamma)
{
    init();
}

RaoControlSurface::~RaoControlSurface()
{
    //dtor
}

void RaoControlSurface::init()
{

    // Database structure
    struct db
    {
        double eta;
        double tht;
        double M;
        double I1;
        double I2;
    };

    // Vector of database points
    std::vector<db> vdb;

    // Angle of the wall at the nozzle exit
    thtE = 0.5*asin(2.0/gamma*(1.0 - Par/isenFlow.Pr(ME))*sqrt(ME*ME-1.0)/(ME*ME));

    // Mach angle at exit of the wall of the nozzle
    muE = MoC.mu(ME);

    // CE - rhs of eq.(17) of Rao's paper
    double CE = cos(thtE-muE)/cos(muE)*1.0/sqrt((gamma-1.0)+2.0/(ME*ME));

    // Initializing min and max values (to be calculated later)
    Mmin   = 1E10;
    Mmax   = 0.0;
    thtmin = 1E10;
    thtmax = 0.0;

    // Adding the first point of the database
    {
        db p;

        p.eta =  1.0;
        p.tht = thtE;
        p.M   =   ME;
        p.I1  =  0.0;
        p.I2  =  0.0;

        if ( p.tht  < thtmin ) thtmin = p.tht;
        if ( thtmax < p.tht  ) thtmax = p.tht;
        if ( p.M    < Mmin   ) Mmin   = p.M;
        if ( Mmax   < p.M    ) Mmax   = p.M;

        //printf("%f %f %f %f %f\n",p.tht*180/acos(-1.0), p.eta, p.M, p.I1, p.I2);

        vdb.push_back(p);
    }

    /*
        Dividing the interval of tht
    */
    {
        double dtht = thtE / ((double) N);

        double tht = thtE;

        while ( tht < acos(-1.0)/2 )
        {
            db p;

            auto& p0 = vdb.back();

            tht += dtht;

            p.tht = tht;

            p.M   = Mt(p.tht);

            p.eta = eta(p.tht, p.M);

            double deta = p.eta - p0.eta;

            if ( deta >= 0.0 ) break;

            if ( tht    < thtmin ) thtmin = tht;
            if ( thtmax < tht    ) thtmax = tht;
            if ( p.M    < Mmin   ) Mmin   = p.M;
            if ( Mmax   < p.M    ) Mmax   = p.M;

            double I1a = isenFlow.rhor(p0.M) * sqrt(isenFlow.Tr(p0.M)) / sin(p0.tht + MoC.mu(p0.M)) * p0.eta;
            double I1b = isenFlow.rhor(p.M ) * sqrt(isenFlow.Tr(p.M )) / sin(p.tht  + MoC.mu(p.M )) * p.eta;

            p.I1 = p0.I1 - deta * (I1a+I1b) / 2;

            double I2a = 1.0/tan(p0.tht + MoC.mu(p0.M));
            double I2b = 1.0/tan(p.tht  + MoC.mu(p.M ));

            p.I2 = p0.I2 - deta * (I2a+I2b) / 2;

            //printf("%f %f %f %f %f\n", p.tht*180/acos(-1.0), p.eta, p.M, p.I1, p.I2);

            vdb.push_back(p);
        }

    }

    // Sorting the database vector from lowest to greatest eta
    std::sort(vdb.begin(),vdb.end(),[](const db& a, const db& b)
    {
        return a.eta < b.eta;
    });

    // Finding minimum value of eta
    etamin = vdb.front().eta;

    std::vector<double> veta;
    std::vector<double> vtht;
    std::vector<double> vI1;
    std::vector<double> vI2;

    for (auto p: vdb)
    {
        veta.push_back(p.eta);
        vtht.push_back(p.tht);
        vI1.push_back(p.I1);
        vI2.push_back(p.I2);

        // Printing data
        //printf("%f %f %f %f %f\n", p.tht, p.eta, p.M, p.I1, p.I2);
        //printf("%f %f %f %f %f\n", p.tht*180/acos(-1.0), p.eta, p.M, p.I1, p.I2);
    }
    //printf("%d\n",vdb.size());

    // Initializing interpolators
    thtInterp.init(veta,vtht);

    I1Interp.init(veta,vI1);

    I2Interp.init(veta,vI2);

    // Minimum value of M
    {
        double raux = -(CE*sqrt(CE*CE*gamma*gamma+2.0*CE*CE*gamma+CE*CE-8.0))
                      /(CE*CE*gamma-CE*CE-1.0)+(CE*CE*gamma)/(CE*CE*gamma-CE*CE-1.0)
                      -(3.0*CE*CE)/(CE*CE*gamma-CE*CE-1.0);

        Mmin = (raux > 0 ? sqrt(raux)/sqrt(2.0) : Mmin);
    }


//    printf("thtmin %f\n", thtmin);
//    printf("thtmax %f\n", thtmax);
//    printf("Mmin   %f\n", Mmin);
//    printf("Mmax   %f\n", Mmax);
}


double RaoControlSurface::Mt(const double & tht, size_t Nit, double tol)
{

    double M0 = ME;
    double M  = M0;

    for (size_t i = 1; i < Nit; ++i)
    {
        double F = cos(thtE-muE)/(sqrt(gamma+2.0/(ME*ME)-1.0)*cos(muE))
                   -cos(tht-asin(1.0/M))/(sqrt(1.0-1.0/(M*M))*sqrt(gamma+2.0/(M*M)-1.0));

        double dFdM = sin(tht-asin(1.0/M))
                      /((1.0-1.0/(M*M))*(M*M)*sqrt(gamma+2.0/(M*M)-1.0))
                      + cos(tht-asin(1.0/M))/(pow(1-1/(M*M),1.5)*M*M*M*sqrt(gamma+2.0/(M*M)-1))
                      -(2.0*cos(tht-asin(1.0/M)))/(sqrt(1.0-1.0/(M*M))*M*M*M*pow(gamma+2.0/(M*M)-1.0,1.5));

        M = M - F/dFdM;

        if ( std::abs(M-M0) < tol ) return M;

        M0 = M;
    }

    printf("RaoControlSurface::Mt: convergence failure!\n");

    return ME;
}


double RaoControlSurface::eta(const double & tht, const double & M)
{
    return g(thtE, ME)/g(tht,M);
}



double RaoControlSurface::g(const double & tht, const double & M)
{
    return M*M / sqrt(M*M-1.0) * isenFlow.Pr(M) * pow(sin(tht), 2.0);
}


double RaoControlSurface::residual(const std::vector<MoCPoint>& cminus)
{
    MoCPoint pointD;

    for (size_t i = 1; i < cminus.size(); ++i)
    {

        /*

            Given the two successive points p1 and p2 along the C-. In this interval
            the flow angle 'theta' varies between theta1 and theta2. Does this interval
            intersect the interval of theta [thtmin, thtmax] of the control surface?
            Lets check it...

        */
        double tht1 = cminus[i-1].tht;
        double tht2 = cminus[i  ].tht;

        double tlower; // Lower value of theta common to C- and to the control surface
        double tupper; // Upper value of theta common to C- and to the control surface

        if ( ! isThereACommonIntervalOfTheta(tht1, tht2, tlower, tupper) ) continue;

        /*

            The C- between p1 and p2 is interpolated by a straight line.
            To find the intersection of this line with the control surface,
            M(tht) of C- must be equal Mt(tht) of the control surface.
            Given func(tht) = M(tht)-Mt(tht), a solution exists if

                        func(tlower) * func(tupper) <= 0

            Lets check it...
        */

        // M and tht of two successive points along the C-
        double M1   = cminus[i-1].M;
        double M2   = cminus[i  ].M;

        // Lambda function to find thtD
        auto func = [&](const double& tht)
        {
            return (M2-M1) * (tht-tht1) / (tht2-tht1) + M1 - Mt(tht);
        };

        /*
            Is there a root within [tlower, tupper]?
        */
        if ( func(tlower) * func(tupper) > 0.0 ) continue;

        // Finding tht of point D
        auto thtD = NumRootFinding::bisection(func, tlower, tupper);

        // Interpolating the remaining quantities at point D
        auto& p1 = cminus[i-1];
        auto& p2 = cminus[i];

        pointD = MoC.linear_interp_tht(p1, p2, thtD);

        // Creating the C- line from B to D
        CmBD.clear();

        for (size_t j = 0; j < i; ++j)
        {
            CmBD.push_back(cminus[j]);
        }

        CmBD.push_back(pointD);

        double etaD = eta(pointD.tht, pointD.M);

        if ( etaD < etamin ) {

            printf("RaoControlSurface::residual error: eta of point D must be greater than eta_min\n");
            return 1E10;
        }

        rE = pointD.r / etaD;

        /*
            Preparing to calculate the integral I3 along BD
        */
        std::vector<double> vx (CmBD.size());
        std::vector<double> vf3(CmBD.size());

        for (size_t j = 0; j < CmBD.size(); ++j)
        {
            double x   = CmBD[j].x;
            double r   = CmBD[j].r;
            double M   = CmBD[j].M;
            double rho = isenFlow.rhor(M);
            double T   = isenFlow.Tr(M);
            double tht = CmBD[j].tht;
            double mu  = CmBD[j].mu;

            vx[j] = x;

            vf3[j] = rho * sqrt(T) * r / cos(tht-mu);
        }

        double I3 = NumIntegration::int_trapezoidal_rule(vf3, vx);

        double res = -I3 + rE * rE * I1Interp.eval(etaD);

        xE = pointD.x + rE * I2Interp.eval(etaD);

//        printf("     res: %f\n",res);
//        printf("      I3: %f\n",I3);
//        printf("rE**2*I1: %f\n",rE * rE * I1Interp.eval(etaD));
//        printf("    thtD: %f\n",pointD.tht);
//        printf("      MD: %f\n",pointD.M);
//        printf("      xD: %f\n",pointD.x);
//        printf("      rD: %f\n",pointD.r);
//        printf("      xE: %f\n",xE);
//        printf("      rE: %f\n",rE);

        return res;
    }

    return 1000.0;
}

std::vector<MoCPoint> RaoControlSurface::cplusDE(size_t N)
{

    // Point D
    MoCPoint pointD = CmBD.back();

    // C+ from point D to point E
    std::vector<MoCPoint> CpDE;

    // Adding point D to CpDE
    CpDE.push_back(pointD);

    // eta at point D
    double etaD = pointD.r / rE;

    // Number of partitions
    size_t Np = N-1;

    // Partitioning the interval of eta
    double deta = (1.0-etaD)/((double) Np);

    // Jumping first point, since it was already added (point D)
    for (size_t i = 1; i <= Np; ++i)
    {
        MoCPoint p;

        double eta = etaD + deta * i;

        p.x = xE - rE * I2Interp.eval(eta);

        p.r = rE * eta;

        p.tht = thtInterp.eval(eta);

        p.M = Mt(p.tht);

        p.mu = MoC.mu(p.M);

        p.nu = MoC.nu(p.M);

        CpDE.push_back(std::move(p));
    }

    return CpDE;
}
