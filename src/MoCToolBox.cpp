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

#include "MoCToolBox.h"
#include <cmath>
#include <cstdio>

MoCToolBox::MoCToolBox(double gamma, int Nit, double eps, double mzero):
    gamma{gamma}, Nit{Nit}, eps{eps}, mzero{mzero}
{
    //ctor
}

MoCToolBox::~MoCToolBox()
{
    //dtor
}


double MoCToolBox::nu(const double& M)
{
    return sqrt((gamma+1.0)/(gamma-1.0))*atan(sqrt((gamma-1.0)/(gamma+1.0)*(M*M-1.0)))-atan(sqrt(M*M-1.0));
}


double MoCToolBox::nu_inv(const double& nu_val)
{
    double Mi {1.1};
    double Mf {1.1};

    for (int i = 0; i < 1000; ++i)
    {

        Mf = Mi - (nu(Mi)-nu_val)/dnudM(Mi);

        //printf("%d  %14.7f\n", i, Mf);

        if ( std::abs(Mf-Mi) < eps ) break;

        Mi = Mf;

    }

    return Mf;
}


double MoCToolBox::dnudM(const double& M)
{
    return (M*(gamma-1)*sqrt((gamma+1)/(gamma-1)))/((gamma+1)*sqrt(((M*M-1)*(gamma-1))/(gamma+1))*(((M*M-1)*(gamma-1))/(gamma+1)+1))-1/(M*sqrt(M*M-1));
}


double MoCToolBox::mu(const double& M)
{
    return asin(1.0/M);
}


double MoCToolBox::fp(const double& tht,const  double& mu)
{
    return tan(tht+mu);
}


double MoCToolBox::fm(const double& tht,const  double& mu)
{
    return tan(tht-mu);
}


double MoCToolBox::gp(const double& tht, const double& M, const double& r)
{
    return ( r > mzero ? -1.0/(r*(sqrt(M*M-1.0)+1.0/tan(tht))): 0.0);
}


double MoCToolBox::gm(const double& tht, const double& M, const double& r)
{
    return ( r > mzero ? 1.0/(r*(sqrt(M*M-1.0)-1.0/tan(tht))) : 0.0 );
}


/*
    p1 is a point lying on the C- characteristics
    p2 is a point lying on the C+ characteristics
*/
MoCPoint MoCToolBox::Cminus_to_Cplus(const MoCPoint& p1, const MoCPoint& p2)
{

    MoCPoint p3;

    // Auxiliary variables
    double tht13;
    double  mu13;
    double   M13;
    double   r13;
    double tht23;
    double  mu23;
    double   M23;
    double   r23;
    double  fm13;
    double  fp23;
    double  gm13;
    double  gp23;

    // Initializing p3
    p3.x   = 0;
    p3.r   = 0;
    p3.tht = 0;
    p3.nu  = 0;
    p3.M   = 0;
    p3.mu  = 0;

    // First order approximation
    tht13 =  p1.tht;
    mu13  =  p1.mu;
    M13   =  p1.M;
    r13   =  p1.r;
    tht23 =  p2.tht;
    mu23  =  p2.mu;
    M23   =  p2.M;
    r23   =  p2.r;

    fm13 = fm(tht13, mu13);

    fp23 = fp(tht23, mu23);

    p3.x   = ( p2.r - p1.r + p1.x * fm13 - p2.x * fp23 )/( fm13 - fp23 );

    p3.r   = p1.r + fm13 * (p3.x-p1.x);

    // If the point p1 has r=0
    if ( std::abs(r13) < mzero )
    {

        gp23 = gp(tht23, M23, r23);

        p3.tht = 1.0/3.0*(p1.nu+p2.tht-p2.nu+gp23*(p3.r-p2.r));

        p3.nu  = -p2.tht+p2.nu+p3.tht-gp23*(p3.r-p2.r);

        // If the point p2 has r=0
    }
    else if ( std::abs(r23) < mzero )
    {

        gm13 = gm(tht13, M13, r13);

        p3.tht = 1.0/3.0*(p1.tht+p1.nu-p2.nu+gm13*(p3.r-p1.r));

        p3.nu  = p1.tht+p1.nu-p3.tht+gm13*(p3.r-p1.r);

        // Otherwise
    }
    else
    {

        gm13 = gm(tht13, M13, r13);

        gp23 = gp(tht23, M23, r23);

        p3.tht = 0.5*(p1.tht+p1.nu+p2.tht-p2.nu+gm13*(p3.r-p1.r)+gp23*(p3.r-p2.r));

        p3.nu  = 0.5*(p1.tht+p1.nu-p2.tht+p2.nu+gm13*(p3.r-p1.r)-gp23*(p3.r-p2.r));
    }

    p3.M   = nu_inv(p3.nu);

    p3.mu  = mu(p3.M);


    //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p3.x, p3.r, p3.tht, p3.nu, p3.M, p3.mu);

    double M0 = p3.M;

    // Second order approximation
    for (int i = 0; i < Nit; ++i)
    {

        tht13 = (p1.tht + p3.tht)/2;
        mu13  = (p1.mu  + p3.mu )/2;
        M13   = (p1.M   + p3.M  )/2;
        r13   = (p1.r   + p3.r  )/2;
        tht23 = (p2.tht + p3.tht)/2;
        mu23  = (p2.mu  + p3.mu )/2;
        M23   = (p2.M   + p3.M  )/2;
        r23   = (p2.r   + p3.r  )/2;

        fm13 = fm(tht13, mu13);
        fp23 = fp(tht23, mu23);
        gm13 = gm(tht13, M13, r13);
        gp23 = gp(tht23, M23, r23);

        p3.x   = ( p2.r - p1.r + p1.x * fm13 - p2.x * fp23 )/( fm13 - fp23 );

        p3.r   = p1.r + fm13 * (p3.x-p1.x);

        p3.tht = 0.5*(p1.tht+p1.nu+p2.tht-p2.nu+gm13*(p3.r-p1.r)+gp23*(p3.r-p2.r));

        p3.nu  = 0.5*(p1.tht+p1.nu-p2.tht+p2.nu+gm13*(p3.r-p1.r)-gp23*(p3.r-p2.r));

        p3.M   = nu_inv(p3.nu);

        p3.mu  = mu(p3.M);

        if ( std::abs(p3.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("Cminus_to_Cplus: convergence failure... Residual of M: %23.16e\n",std::abs(p3.M-M0));
//            printp(p1);
//            printp(p2);
//            printp(p3);
        }

        M0 = p3.M;

    }

    return std::move(p3);
}


// p1 is a point lying on the C+ characteristic
MoCPoint MoCToolBox::Cplus_to_wall(Wall& wall, const MoCPoint& p1)
{

    MoCPoint p2;

    // Auxiliary variables
    double     a;
    double     b;
    double tht12;
    double  mu12;
    double   M12;
    double   r12;
    double  fp12;
    double  gp12;

    // Initializing p2
    p2 = p1;

    p2.r = wall.r(p2.x);

    double M0 = p2.M;

    // Second order approximation
    for (int i = 0; i < Nit; ++i)
    {

        wall.linear_approximation(p2.x, a, b);

        p2.tht = atan(a);

        tht12 = (p1.tht + p2.tht)/2;
        mu12  = (p1.mu  + p2.mu )/2;
        M12   = (p1.M   + p2.M  )/2;
        r12   = (p1.r   + p2.r  )/2;

        fp12 = fp(tht12, mu12);
        gp12 = gp(tht12, M12, r12);

        p2.x   = (p1.r-fp12*p1.x-b)/(a-fp12);

        p2.r   = a*p2.x+b;

        p2.nu  = p2.tht-p1.tht+p1.nu-gp12*(p2.r-p1.r);

        p2.M   = nu_inv(p2.nu);

        p2.mu  = mu(p2.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);


        if ( std::abs(p2.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("\n Cplus_to_wall: convergence failure...\n");
            printf("Residual of M: %g \n",std::abs(p2.M-M0));
        }

        M0 = p2.M;

    }

    return std::move(p2);
}


// p1 is a point lying on the C- characteristic
MoCPoint MoCToolBox::Cminus_to_wall(Wall& wall, const MoCPoint& p1)
{

    MoCPoint p2;

    // Auxiliary variables
    double     a;
    double     b;
    double tht12;
    double  mu12;
    double   M12;
    double   r12;
    double  fm12;
    double  gm12;

    // Initializing p2
    p2   = p1;
    p2.r = wall.r(p2.x);

    double M0 = p2.M;

    // Second order approximation
    for (int i = 0; i < Nit; ++i)
    {

        wall.linear_approximation(p2.x, a, b);

        p2.tht = atan(a);

        tht12 = (p1.tht + p2.tht)/2;
        mu12  = (p1.mu  + p2.mu )/2;
        M12   = (p1.M   + p2.M  )/2;
        r12   = (p1.r   + p2.r  )/2;

        fm12 = fm(tht12, mu12);
        gm12 = gm(tht12, M12, r12);

        p2.x   = (p1.r-fm12*p1.x-b)/(a-fm12);

        p2.r   = a*p2.x+b;

        p2.nu  = p1.tht+p1.nu-p2.tht+gm12*(p2.r-p1.r);

        p2.M   = nu_inv(p2.nu);

        p2.mu  = mu(p2.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);


        if ( std::abs(p2.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("\n Cminus_to_wall: convergence failure...\n");
            printf("Residual of M: %g \n",std::abs(p2.M-M0));
        }

        M0 = p2.M;
    }

    return std::move(p2);
}


// p1 is a point lying on the C- characteristic
MoCPoint MoCToolBox::Cminus_to_symmetry_line(const MoCPoint& p1)
{

    MoCPoint p2;

    // Auxiliary variables
    double tht12;
    double  mu12;
    double   M12;
    double   r12;
    double  fm12;
    double  gm12;

    // Initializing p2
    p2.x   = 0;
    p2.r   = 0;
    p2.tht = 0;
    p2.nu  = 0;
    p2.M   = 0;
    p2.mu  = 0;

    // First order approximation
    tht12 =  p1.tht;
    mu12  =  p1.mu;
    M12   =  p1.M;
    r12   =  p1.r;

    fm12 = fm(tht12, mu12);
    gm12 = gm(tht12, M12, r12);

    p2.x   = p1.x - p1.r/fm12;

    p2.nu  = p1.tht + p1.nu - gm12 * p1.r;

    p2.M   = nu_inv(p2.nu);

    p2.mu  = mu(p2.M);

    //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);

    double M0 = p2.M;

    // Second order approximation
    for (int i = 0; i < Nit; ++i)
    {

        tht12 = (p1.tht + p2.tht)/2;
        mu12  = (p1.mu  + p2.mu )/2;
        M12   = (p1.M   + p2.M  )/2;
        r12   = (p1.r   + p2.r  )/2;

        fm12 = fm(tht12, mu12);
        gm12 = gm(tht12, M12, r12);

        p2.x   = p1.x - p1.r/fm12;

        p2.nu  = p1.tht + p1.nu - gm12 * p1.r;

        p2.M   = nu_inv(p2.nu);

        p2.mu  = mu(p2.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);

        if ( std::abs(p2.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("Cminus_to_symmetry_line: convergence failure... Residual of M: %g \n",std::abs(p2.M-M0));
        }

        M0 = p2.M;
    }

    return std::move(p2);
}


// p1 is a point lying on the C+ characteristic
MoCPoint MoCToolBox::Cplus_to_symmetry_line(const MoCPoint& p1)
{

    MoCPoint p2;

    // Auxiliary variables
    double tht12;
    double  mu12;
    double   M12;
    double   r12;
    double  fp12;
    double  gp12;

    // Initializing p2
    p2.x   = 0;
    p2.r   = 0;
    p2.tht = 0;
    p2.nu  = 0;
    p2.M   = 0;
    p2.mu  = 0;

    // First order approximation
    tht12 =  p1.tht;
    mu12  =  p1.mu;
    M12   =  p1.M;
    r12   =  p1.r;

    fp12 = fp(tht12, mu12);
    gp12 = gp(tht12, M12, r12);

    p2.x   = p1.x - p1.r/fp12;

    p2.nu  = gp12 * p1.r + p1.nu - p1.tht;

    p2.M   = nu_inv(p2.nu);

    p2.mu  = mu(p2.M);

    //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);

    double M0 = p2.M;

    // Second order approximation
    for (int i = 0; i < Nit; ++i)
    {

        tht12 = (p1.tht + p2.tht)/2;
        mu12  = (p1.mu  + p2.mu )/2;
        M12   = (p1.M   + p2.M  )/2;
        r12   = (p1.r   + p2.r  )/2;

        fp12 = fp(tht12, mu12);
        gp12 = gp(tht12, M12, r12);

        p2.x   = p1.x - p1.r/fp12;

        p2.nu  = gp12 * p1.r + p1.nu - p1.tht;

        p2.M   = nu_inv(p2.nu);

        p2.mu  = mu(p2.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", p2.x, p2.r, p2.tht, p2.nu, p2.M, p2.mu);

        if ( std::abs(p2.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("\n Cplus_to_symmetry_line: convergence failure...\n");
            printf("Residual of M: %g \n",std::abs(p2.M-M0));
        }

        M0 = p2.M;
    }

    return std::move(p2);
}


MoCPoint MoCToolBox::wall_to_Cminus(const MoCPoint& p1, const MoCPoint& p2, MoCPoint& pw)
{

    MoCPoint p3;

    double zeta;
    double r3w;
    double M3w;
    double tht3w;
    double mu3w;
    double fp3w;
    double gp3w;

    /*

        First approximation

    */
    zeta   = 0.5;

    p3.x   = (p2.x  -p1.x  ) * zeta + p1.x;
    p3.r   = (p2.r  -p1.r  ) * zeta + p1.r;
    p3.M   = (p2.M  -p1.M  ) * zeta + p1.M;
    p3.tht = (p2.tht-p1.tht) * zeta + p1.tht;
    p3.mu  = mu(p3.M);
    p3.nu  = nu(p3.M);

    /*

        pw.x   --> Given
        pw.r   --> Given
        pw.tht --> Given

    */
    pw.M   = p3.M;
    pw.mu  = p3.mu;
    pw.nu  = p3.nu;

    double M0 = pw.M;

    for (int i = 0; i < Nit; ++i)
    {

        tht3w = ( p3.tht + pw.tht ) / 2.0;
        mu3w  = ( p3.mu  + pw.mu  ) / 2.0;

        fp3w = fp(tht3w, mu3w);

        zeta = (fp3w * (pw.x-p1.x) - (pw.r-p1.r)) / (fp3w * (p2.x-p1.x) - (p2.r-p1.r));

        p3.x   = (p2.x  -p1.x  ) * zeta + p1.x;
        p3.r   = (p2.r  -p1.r  ) * zeta + p1.r;
        p3.M   = (p2.M  -p1.M  ) * zeta + p1.M;
        p3.tht = (p2.tht-p1.tht) * zeta + p1.tht;
        p3.mu  = mu(p3.M);
        p3.nu  = nu(p3.M);

        tht3w = ( p3.tht + pw.tht ) / 2.0;
        M3w   = ( p3.M   + pw.M   ) / 2.0;
        r3w   = ( p3.r   + pw.r   ) / 2.0;

        gp3w = gp(tht3w, M3w, r3w);

        pw.nu = pw.tht - p3.tht + p3.nu - gp3w * ( pw.r - p3.r );

        pw.M = nu_inv(pw.nu);

        pw.mu = mu(pw.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", pw.x, pw.r, pw.tht, pw.nu, pw.M, pw.mu);
        //printf("%g\n",zeta);

        if ( std::abs(pw.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("\n wall_to_Cminus: convergence failure...\n");
            printf("Residual of M: %g \n",std::abs(pw.M-M0));
        }

        M0 = pw.M;

    }

    return std::move(p3);
}


MoCPoint MoCToolBox::wall_to_Cplus(const MoCPoint& p1, const MoCPoint& p2, MoCPoint& pw)
{

    MoCPoint p3;

    double zeta;
    double r3w;
    double M3w;
    double tht3w;
    double mu3w;
    double fm3w;
    double gm3w;

    /*

        First approximation

    */
    zeta   = 0.5;

    p3.x   = (p2.x  -p1.x  ) * zeta + p1.x;
    p3.r   = (p2.r  -p1.r  ) * zeta + p1.r;
    p3.M   = (p2.M  -p1.M  ) * zeta + p1.M;
    p3.tht = (p2.tht-p1.tht) * zeta + p1.tht;
    p3.mu  = mu(p3.M);
    p3.nu  = nu(p3.M);

    /*

        pw.x   --> Given
        pw.r   --> Given
        pw.tht --> Given

    */
    pw.M   = p3.M;
    pw.mu  = p3.mu;
    pw.nu  = p3.nu;

    double M0 = pw.M;

    for (int i = 0; i < Nit; ++i)
    {

        tht3w = ( p3.tht + pw.tht ) / 2.0;
        mu3w  = ( p3.mu  + pw.mu  ) / 2.0;

        fm3w = fm(tht3w, mu3w);

        zeta = (fm3w * (pw.x-p1.x) - (pw.r-p1.r)) / (fm3w * (p2.x-p1.x) - (p2.r-p1.r));

        p3.x   = (p2.x  -p1.x  ) * zeta + p1.x;
        p3.r   = (p2.r  -p1.r  ) * zeta + p1.r;
        p3.M   = (p2.M  -p1.M  ) * zeta + p1.M;
        p3.tht = (p2.tht-p1.tht) * zeta + p1.tht;
        p3.mu  = mu(p3.M);
        p3.nu  = nu(p3.M);

        tht3w = ( p3.tht + pw.tht ) / 2.0;
        M3w   = ( p3.M   + pw.M   ) / 2.0;
        r3w   = ( p3.r   + pw.r   ) / 2.0;

        gm3w = gm(tht3w, M3w, r3w);

        pw.nu = -pw.tht + p3.tht + p3.nu + gm3w * ( pw.r - p3.r );

        pw.M = nu_inv(pw.nu);

        pw.mu = mu(pw.M);

        //printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", pw.x, pw.r, pw.tht, pw.nu, pw.M, pw.mu);

        if ( std::abs(pw.M-M0) < eps ) break;

        if ( i == Nit-1 )
        {
            printf("\n wall_to_Cplus: convergence failure...\n");
            printf("Residual of M: %g \n",std::abs(pw.M-M0));
        }

        M0 = pw.M;

    }

    return std::move(p3);
}


void printp(const MoCPoint& p, bool header)
{
    if (header) printf("#%9s, %10s, %10s, %10s, %10s, %10s\n","x","r","tht","M","nu(M)","mu(M)");
    printf("%10.5g, %10.5g, %10.5g, %10.5g, %10.5g, %10.5g\n", p.x, p.r, p.tht, p.M, p.nu, p.mu);
}
