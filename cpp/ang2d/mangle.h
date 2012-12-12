/*
 * mangle.h
 *
 *   Extracted from the make_random.cpp file by Martin White.
 *
 *  Created on: Dec 12, 2012
 *      Author: npadmana
 */

#ifndef MANGLE_H_
#define MANGLE_H_

#include <vector>
#include <list>
#include <string>

namespace Mangle {

/** Class definition for a cap.
 *
 * See mangle documentation for a detailed description of what all this means.
 */
class CapClass {
private:
	double	x,y,z,cm;
public:
	// Constructor
	CapClass(double x1, double y1, double z1, double cm1) : x(x1), y(y1), z(z1), cm(cm1) {};

	// Test if in cap
	bool incap(double x0, double y0, double z0);
	bool incap(double theta, double phi);
};	// CapClass


/** Class for Mangle polygons
 *
 * See mangle documentation for more.
 */
class PolygonClass {
private:
	std::vector<CapClass>	caps;
	long			pid;
	double		weight;
public:
	PolygonClass() : weight(0.0), pid(0) {};

	// Getters and setters
	long getid();
	void setid(long pid1);
	double getwt();
	void setwt(double wt);

	// Add in a cap
	void addcap(CapClass cap1);

	// Check to see if in polygons using x,y,z
	bool inpoly(double x, double y, double z);

	// Check to see if in polygons using theta, phi
	bool inpoly(double theta, double phi) ;


};	// PolygonClass

class SDSSpixClass {	// Does conversions for us.  Based on SDSSpix.pro
private:
	long		nx0,ny0,resolution;
	double	deg2Rad,node,etapole;
	double	etaOffSet,surveyCenterRA,surveyCenterDEC;
	std::vector<double> Mangle::SDSSpixClass::tp2survey(double theta, double phi);
public:
	SDSSpixClass(long pixres);
	long Mangle::SDSSpixClass::pixelnum(double theta, double phi);
};	// SDSSpixClass

class MaskClass {
// A class to handle masks.  When invoked it reads an ascii mask file.
// The main method returns a completeness, used for making random catalogs.
private:
  std::vector<PolygonClass>	polygons;
  std::vector< std::list<int> >	pixels;
  int				pixelres;
  char				pixeltype;
  long parsepoly(std::string sbuf, long &ncap, double &weight, long &pixel);
  long pixelnum(double theta, double phi);
public:
  MaskClass(std::string fname);
  long npolygons() {	// For general interest, how many polygons in mask
    return(polygons.size());
  }
  double completeness(double theta, double phi) {
  // This is the main method, returning the completness at (theta,phi).
    bool	notfnd=true;
    long	ii,pid=-1;
    double	wt=0;
    if (pixelres==-1) {	// Mask isn't pixelized.
      for (ii=0; ii<polygons.size() && notfnd; ii++)
        if (polygons[ii].inpoly(theta,phi)) {
          wt    = polygons[ii].getwt();
          notfnd= false;
        }
    }
    else {
      long ipix;
      if (pixeltype=='s')
        ipix= this->pixelnum(theta,phi);
      else if (pixeltype=='d') {
        SDSSpixClass sdsspix(pixelres);
        ipix= sdsspix.pixelnum(theta,phi);
      }
      else {
        cerr << "Unknown pixelization scheme " << pixeltype << endl;
        exit(1);
      }
      if (ipix<pixels.size()) {
        for (std::list<int>::iterator ii=pixels[ipix].begin();
             ii!=pixels[ipix].end() && notfnd; ii++) {
          if (polygons[*ii].inpoly(theta,phi)) {
            wt    = polygons[*ii].getwt();
            notfnd= false;
          }
        }
      }
    }
    if (notfnd)
      return(0.0);
    else
      return(wt);
  }
};




} // end namespace mangle




#endif /* MANGLE_H_ */
