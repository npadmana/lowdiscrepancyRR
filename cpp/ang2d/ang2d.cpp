#include "ang2d.h"

void Ang2D::vecshift(vector<double> &x, vector<double> &x0) {
	double tmp;
	transform(x0.begin(), x0.end(), x.begin(), x.begin(),
			[&tmp](double a, double b){return modf(a+b,&tmp);});
}
