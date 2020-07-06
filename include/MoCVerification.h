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

#ifndef MOCVERIFICATION_H
#define MOCVERIFICATION_H

#include "MoCToolBox.h"
#include "Wall.h"
#include <cmath>
#include <cstdio>

/**
    \brief MoCVerification performs calculations for MoCToolBox code verification
*/
class MoCVerification
{

public:

    /**
        \brief Constructor
    */
    MoCVerification():
        moc(1.4)
    {

        p1.x   = 0.1;
        p1.r   = 0.3;
        p1.M   = 1.9;
        p1.tht = 20.0 * acos(-1.0) / 180.0;
        p1.nu  = moc.nu(p1.M);
        p1.mu  = moc.mu(p1.M);

        p2.x   = 0.12;
        p2.r   = 0.08;
        p2.M   = 1.8;
        p2.tht = 10.0 * acos(-1.0) / 180.0;
        p2.nu  = moc.nu(p2.M);
        p2.mu  = moc.mu(p2.M);

        p3.x   = 0.12;
        p3.r   = 0.0;
        p3.M   = 1.8;
        p3.tht = 0.0;
        p3.nu  = moc.nu(p3.M);
        p3.mu  = moc.mu(p3.M);

        pw.x   = 0.12;
        pw.r   = uwall.r(pw.x);
        pw.tht = atan(uwall.drdx(pw.x));
//        pw.M   = ??;
//        pw.nu  = ??;
//        pw.mu  = ??;

    }


    /**
        \brief Calculates 10 verification cases
    */
    void calculate(){
        MoCPoint p;

        printf("\n Test 1: Cminus_to_Cplus(p1,p2)\n");
        p = moc.Cminus_to_Cplus(p1,p2);
        printf("Point p1 of C-, point p2 of C+ and point p at C- and C+ intersection:\n");
        printPoint(p1, true);
        printPoint(p2);
        printPoint(p);


        printf("\n Test 2: Cminus_to_Cplus(p2,p1)\n");
        p = moc.Cminus_to_Cplus(p2,p1);
        printf("Point p2 of C-, point p1 of C+ and point p at C- and C+ intersection:\n");
        printPoint(p2, true);
        printPoint(p1);
        printPoint(p);


        printf("\n Test 3: Cminus_to_Cplus(p1,p3) (p3 is on the symmetry line)\n");
        p = moc.Cminus_to_Cplus(p1,p3);
        printf("Point p1 of C-, point p3 of C+ and point p at C- and C+ intersection:\n");
        printPoint(p1, true);
        printPoint(p3);
        printPoint(p);


        printf("\n Test 4: Cminus_to_Cplus(p3,p1) (p3 is on the symmetry line)\n");
        p = moc.Cminus_to_Cplus(p3,p1);
        printf("Point p3 of C-, point p1 of C+ and point p at C- and C+ intersection:\n");
        printPoint(p3, true);
        printPoint(p1);
        printPoint(p);


        printf("\n Test 5: Cminus_to_symmetry_line(p2)\n");
        p = moc.Cminus_to_symmetry_line(p2);
        printf("Point p2 of C-, point p at symmetry line intersection:\n");
        printPoint(p2, true);
        printPoint(p);


        printf("\n Test 6: Cplus_to_symmetry_line(p2)\n");
        p = moc.Cplus_to_symmetry_line(p2);
        printf("Point p2 of C+, point p at symmetry line intersection:\n");
        printPoint(p2, true);
        printPoint(p);


        printf("\n Test 7: Cplus_to_wall(uwall, p1)\n");
        p = moc.Cplus_to_wall(uwall, p1);
        printf("Point p1 of C+, point p at upper wall intersection:\n");
        printPoint(p1, true);
        printPoint(p);


        printf("\n Test 8: Cminus_to_wall(lwall, p2)\n");
        p = moc.Cminus_to_wall(lwall, p2);
        printf("Point p2 of C-, point p at lower wall intersection:\n");
        printPoint(p2, true);
        printPoint(p);


        printf("\n Test 9: wall_to_Cminus(p1, p2, pw)\n");
        p = moc.wall_to_Cminus(p1, p2, pw);
        printf("Points p1 and p2 of C-, point pw on the wall and p (intersection of C- and C+):\n");
        printPoint(p1, true);
        printPoint(p2);
        printPoint(pw);
        printPoint(p);


        printf("\n Test 10: wall_to_Cplus(p1, p2, pw)\n");
        p = moc.wall_to_Cplus(p1, p2, pw);
        printf("Points p1 and p2 of C+, point pw on the wall and p (intersection of C- and C+):\n");
        printPoint(p1, true);
        printPoint(p2);
        printPoint(pw);
        printPoint(p);
    };

private:

    /**
        \brief Prints MoCPoint
    */
    void printPoint(MoCPoint& p, bool withHeader = false){
        if ( withHeader ) printf("%23s, %23s, %23s, %23s, %23s, %23s \n", "x", "r", "tht", "M", "nu", "mu");
        printf("%23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e \n", p.x, p.r, p.tht, p.M, p.nu, p.mu);
    }


    /**
        \brief Defines an upper wall
    */
    class UpperWall: public Wall
    {
    public:

        virtual void linear_approximation(const double& x, double& a, double& b)
        {
            a=drdx(x);
            b=r(x)-drdx(x)*x;
        };

        virtual double r(const double& x)
        {
            return 0.1 * x*x + 0.4;
        };

        virtual double drdx(const double& x)
        {
            return 0.2*x;
        };
    };


    /**
        \brief Defines a lower wall
    */
    class LowerWall: public Wall
    {
    public:
        virtual void linear_approximation(const double& x, double& a, double& b)
        {
            a=drdx(x);
            b=r(x)-drdx(x)*x;
        };

        virtual double r(const double& x)
        {
            return -0.1 * x*x + 0.05;
        };

        virtual double drdx(const double& x)
        {
            return -0.2*x;
        };
    };

    UpperWall uwall;
    LowerWall lwall;

    MoCPoint p1;
    MoCPoint p2;
    MoCPoint p3;
    MoCPoint pw;

    MoCToolBox moc;
};

#endif // MOCVERIFICATION_H
