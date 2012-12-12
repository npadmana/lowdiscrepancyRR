#include "mangle.h"

#include	<cmath>
#include	<iostream>
#include	<fstream>
#include	<sstream>
#include	<iomanip>
#include	<vector>
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

long Mangle::PolygonClass::getid() {
	return pid;
}

void Mangle::PolygonClass::setid(long int pid1) {
	pid = pid1;
}

void Mangle::PolygonClass::addcap(CapClass cap1) {
	caps.push_back(cap1);
}

double Mangle::PolygonClass::getwt() {
	return weight;
}

void Mangle::PolygonClass::setwt(double wt) {
	weight = wt;
}

// Note that both these methods use short circuits to jump out of the loop
// as soon as a single polygon fails
bool Mangle::PolygonClass::inpoly(double theta, double phi) {
	bool tmp=true;
	for (std::vector<CapClass>::iterator i=caps.begin();
			i!=caps.end() && tmp; i++) {
		tmp = tmp && ( i->incap(theta,phi) );
	}
	return(tmp);
}

bool Mangle::PolygonClass::inpoly(double x, double y, double z) {
  bool tmp=true;
  for (std::vector<CapClass>::iterator i=caps.begin();
       i!=caps.end() && tmp; i++) {
    tmp = tmp && ( i->incap(x,y,z) );
  }
  return(tmp);
}

std::vector<double> Mangle::SDSSpixClass::tp2survey(double theta, double phi) {
	// Converts (theta,phi) to survey coordinates (eta,lambda).
	// On input (theta,phi) are in radians, eta/lambda come out in degrees.
	std::vector<double> retval(2);
	double x,y,z,ra,dec;

	ra = phi;
	dec= M_PI/2.0-theta;
	x  = cos(ra-node)*cos(dec);
	y  = sin(ra-node)*cos(dec);
	z  = sin(dec);

	double lambda = -1.0*asin(x)/deg2Rad;
	double eta    = (atan2(z,y) - etapole)/deg2Rad;

	if (eta <-180.0) eta += 360.0;
	if (eta > 180.0) eta -= 360.0;
	retval[0] = eta;
	retval[1] = lambda;
	return(retval);
}

Mangle::SDSSpixClass::SDSSpixClass(long pixres) {
	nx0=36;
	ny0=13;
	deg2Rad=M_PI/180.;
	etaOffSet=1.25;
	surveyCenterRA=185.0;
	surveyCenterDEC=32.5;
	node=deg2Rad*(surveyCenterRA-90.0);
	etapole=deg2Rad*surveyCenterDEC;
	resolution = pixres;
}


long Mangle::SDSSpixClass::pixelnum(double theta, double phi) {
	long   i,j,nx=nx0*resolution,ny=ny0*resolution;
	std::vector<double> surv = tp2survey(theta,phi);

	double eta= (surv[0]-etaOffSet)*deg2Rad;
	if (eta<0) eta += 2.*M_PI;
	i = long(nx*eta/2/M_PI);

	double lambda = (90.0-surv[1])*deg2Rad;
	if (lambda>=M_PI)
		j = ny-1;
	else
		j = long( ny*( (1.0-cos(lambda))/2.0 ) );
	return(nx*j+i);
}


long Mangle::MaskClass::parsepoly(std::string sbuf, long &ncap, double &weight, long &pixel) {
	// Reads the "polygon" line in a Mangle file, returning #caps, weight & pixel.
	long	i,j,ipoly=-1;
	std::string	ss;
	if ( (i=sbuf.find("polygon")   )!=std::string::npos &&
			(j=sbuf.find("("))!=std::string::npos) {
		ss.assign(sbuf,i+string("polygon").size(),j);
		istringstream(ss) >> ipoly;
	}
	else {
		cerr << "Cannot parse " << sbuf << endl;
		exit(1);
	}
	if ( (i=sbuf.find("(")   )!=std::string::npos &&
			(j=sbuf.find("caps"))!=std::string::npos) {
		ss.assign(sbuf,i+1,j);
		istringstream(ss) >> ncap;
	}
	else {
		cerr << "Cannot parse " << sbuf << endl;
		exit(1);
	}
	if ( (i=sbuf.find("caps,") )!=std::string::npos &&
			(j=sbuf.find("weight"))!=std::string::npos) {
		ss.assign(sbuf,i+string("caps,").size(),j);
		istringstream(ss) >> weight;
	}
	else {
		cerr << "Cannot parse " << sbuf << endl;
		exit(1);
	}
	if ( (i=sbuf.find("weight,"))!=std::string::npos &&
			(j=sbuf.find("pixel")  )!=std::string::npos) {
		ss.assign(sbuf,i+string("weight,").size(),j);
		istringstream(ss) >> pixel;
	}
	else {
		if (pixeltype=='u' || pixelres<0) {
			pixel = 0;
		}
		else {
			cerr << "Cannot parse " << sbuf << endl;
			exit(1);
		}
	}
	return(ipoly);
}

long Mangle::MaskClass::pixelnum(double theta, double phi) {
	// For the "simple" pixelization we're just Cartesian in cos(theta) & phi.
	long ipix=0;
	if (pixelres>0) {
		long   ps=0,p2=1;
		for (int i=0; i<pixelres; i++) {  // Work out # pixels/dim and start pix.
			p2  = p2<<1;
			ps += (p2/2)*(p2/2);
		}
		double cth = cos(theta);
		long   n   = (cth==1.0)?0:long( ceil( (1.0-cth)/2 * p2 )-1 );
		long   m   = long( floor( (phi/2./M_PI)*p2 ) );
		ipix= p2*n+m + ps;
	}
	return(ipix);
}

Mangle::MaskClass::MaskClass(std::string fname) {	// Load the mask from fname.
	std::string	sbuf;
	long	i,j,npoly,maxpix=-1;
	ifstream	fs(fname.c_str());
	if (!fs) {
		cerr << "Unable to open " << fname << endl;
		exit(1);
	}
	// Read in the number of polygons, which should start the file.
	getline(fs,sbuf);
	if (fs.eof()) {cerr<<"Unexpected end-of-file."<<endl;exit(1);}
	istringstream(sbuf) >> npoly;
	polygons.resize(npoly);
	// See if we're pixelized.
	getline(fs,sbuf);
	if ( (i=sbuf.find("pixelization"))!=std::string::npos &&
			(j=sbuf.find("s",i)         )!=std::string::npos) {
		string ss;
		ss.assign(sbuf,i+string("pixelization").size(),j);
		istringstream(ss) >> pixelres;
		pixeltype = 's';
	}
	else if ( (i=sbuf.find("pixelization"))!=std::string::npos &&
			(j=sbuf.find("d",i)         )!=std::string::npos) {
		string ss;
		ss.assign(sbuf,i+string("pixelization").size(),j);
		istringstream(ss) >> pixelres;
		pixeltype = 'd';
		// It turns out that the Mangle "default" pixelization scheme is
		// not quite SDSSpix, but "based on SDSSpix" in a fairly complicated
		// manner.  So for now, pretend it's not pixelized.
		pixelres  = -1;
		pixeltype = 'u';
	}
	else  {
		pixelres  = -1;
		pixeltype = 'u';
	}
	cout << "# Pixel res "<< pixelres << ", type " << pixeltype << endl;
	cout.flush();
	// For each polygon, create the appropriate "polygons" entry.
	for (int ipoly=0; ipoly<npoly; ipoly++) {
		long	ncap,pixel;
		double	weight;
		do {
			getline(fs,sbuf);
			if (fs.eof()) {cerr<<"Unexpected end-of-file."<<endl;exit(1);}
		} while(sbuf.find("polygon")==std::string::npos);
		j = parsepoly(sbuf,ncap,weight,pixel);
		if (j!=ipoly) {
			cerr<<"Error reading "<<fname<<endl;
			cerr<<"Read polygon "<<j<<" expecting "<<ipoly<<endl;
			exit(1);
		}
		polygons[ipoly].setid(pixel);
		polygons[ipoly].setwt(weight);
		if (pixel>maxpix) maxpix=pixel;
		// and read in the caps.
		for (int icap=0; icap<ncap; icap++) {
			double x1,y1,z1,cm1;
			getline(fs,sbuf);
			if (fs.eof()) {cerr<<"Unexpected end-of-file."<<endl;exit(1);}
			istringstream(sbuf) >> x1 >> y1 >> z1 >> cm1;
			polygons[ipoly].addcap( Mangle::CapClass(x1,y1,z1,cm1) );
		}
	}
	fs.close();
	// If we're pixelized make a list of which polygons lie in each pixel.
	if (pixelres>=0) {
		pixels.resize(maxpix+1);
		for (int ipoly=0; ipoly<npoly; ipoly++)
			pixels[polygons[ipoly].getid()].push_back(ipoly);
	}
}

long Mangle::MaskClass::npolygons() {	// For general interest, how many polygons in mask
	return(polygons.size());
}

double Mangle::MaskClass::completeness(double theta, double phi) {
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

double Mangle::MaskClass::completeness_radec(double ra, double dec) {
	const double d2r = M_PI/180.0;
	double theta, phi;
	phi = ra*d2r;
	theta = (90.0-dec)*d2r;
	return completeness(theta, phi);
}
