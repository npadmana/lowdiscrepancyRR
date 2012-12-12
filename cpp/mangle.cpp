#include "mangle.h"

#include	<cmath>
#include	<iostream>
#include	<fstream>
#include	<sstream>
#include	<iomanip>
#include	<vector>
#include	<list>
#include	<string>

using namespace std;

bool Mangle::CapClass::incap(double x0, double y0, double z0) {
    double	cdot;
    bool	tmp;
    cdot = 1.0 - x*x0 - y*y0 - z*z0;
    if (cm < 0.0)
      tmp = cdot > (-cm);
    else
      tmp = cdot < cm;
    return(tmp);
}

bool Mangle::CapClass::incap(double theta, double phi) {
	double x0, y0, z0;
    x0 = sin(theta)*cos(phi);
    y0 = sin(theta)*sin(phi);
    z0 = cos(theta);

    return incap(x0,y0,z0);
}
