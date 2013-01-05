#include "mangle2.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

bool Mangle2::MaskClass::inCap(const CapClass& cap, const Vector4d &x1) {
    return x1.dot(cap.x) < cap.cm;
    //cdot = 1.0 - cap.x*x0 - cap.y*y0 - cap.z*z0;
    //return (cap.cm < 0.0) ? (cdot > (-cap.cm)) : (cdot < cap.cm);
}

bool Mangle2::MaskClass::inPolygon(const PolygonClass& poly,
		const Vector4d &x1) {
	return all_of(caps.begin()+poly.capbegin,
				  caps.begin()+poly.capbegin+poly.ncaps,
				  [&](const CapClass& c){return inCap(c,x1);});
}






long Mangle2::MaskClass::pixelnum(double theta, double phi) {
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


long Mangle2::MaskClass::npolygons() {	// For general interest, how many polygons in mask
	return(polys.size());
}

int Mangle2::MaskClass::getPixelres() const {
	return pixelres;
}

double Mangle2::MaskClass::findpoly(PolyIterator start, PolyIterator end,
		const Vector4d &x1,
		long & polyid) {
	polyid = -1;
	PolyIterator res = find_if(start, end,
			           [&](const PolygonClass &p)
			           {return inPolygon(p, x1);});
	if (res == end) {
		return 0;
	} else {
		polyid = res->polyid;
		return res->wt;
	}
}

double Mangle2::MaskClass::completeness(double theta, double phi, long &polyid) {
	// This is the main method, returning the completness at (theta,phi).
//	double x0,y0,z0;
//    x0 = sin(theta)*cos(phi);
//    y0 = sin(theta)*sin(phi);
//    z0 = cos(theta);
	Vector4d x1; x1 << 1.0, -sin(theta)*cos(phi), -sin(theta)*sin(phi), -cos(theta);

    if (pixelres==-1) {
    	return findpoly(polys.begin(), polys.end(), x1, polyid);
    } else {
    	long ipix;
    	if (pixeltype=='s') {
    		ipix= this->pixelnum(theta,phi);
    	} else {
    		throw invalid_argument("Unknown pixelization");
    	}

    	// Check to see if we might have an index to this pixel
    	if ((ipix < pixelIndex.size()) && (pixelIndex[ipix] >= 0))  // Short-circuit
    	{
    		auto start = polys.begin()+ pixelIndex[ipix];
    		auto end = start + numPolyInPixel[ipix];
    		return findpoly(start, end, x1, polyid);
    	} else {
    		polyid = -1;
    		return 0.0;
    	}

    }
}



double Mangle2::MaskClass::completeness(double theta, double phi) {
	long polyid;
	return completeness(theta, phi, polyid);
}

double Mangle2::MaskClass::completeness_radec(double ra, double dec, long &polyid) {
	const double d2r = M_PI/180.0;
	double theta, phi;
	phi = ra*d2r;
	theta = (90.0-dec)*d2r;
	return completeness(theta, phi, polyid);
}

double Mangle2::MaskClass::completeness_radec(double ra, double dec) {
	long polyid;
	return completeness_radec(ra,dec, polyid);
}




void Mangle2::MaskClass::setPixelres(int pixelres) {
	this->pixelres = pixelres;
}

long Mangle2::MaskClass::parsepoly(std::string sbuf, long &ncap, double &weight, long &pixel) {
	// Reads the "polygon" line in a Mangle file, returning #caps, weight & pixel.
	long	i,j,ipoly=-1;
	std::string	ss;
	if ( (i=sbuf.find("polygon")   )!=std::string::npos &&
			(j=sbuf.find("("))!=std::string::npos) {
		ss.assign(sbuf,i+string("polygon").size(),j);
		istringstream(ss) >> ipoly;
	}
	else {
		throw runtime_error(string("Cannot parse ")+=sbuf);
	}
	if ( (i=sbuf.find("(")   )!=std::string::npos &&
			(j=sbuf.find("caps"))!=std::string::npos) {
		ss.assign(sbuf,i+1,j);
		istringstream(ss) >> ncap;
	}
	else {
		throw runtime_error(string("Cannot parse ")+=sbuf);
	}
	if ( (i=sbuf.find("caps,") )!=std::string::npos &&
			(j=sbuf.find("weight"))!=std::string::npos) {
		ss.assign(sbuf,i+string("caps,").size(),j);
		istringstream(ss) >> weight;
	}
	else {
		throw runtime_error(string("Cannot parse ")+=sbuf);
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
			throw runtime_error(string("Cannot parse ")+=sbuf);
		}
	}
	return(ipoly);
}

Mangle2::MaskClass::MaskClass(string fname) :
				pixeltype('u'),
				pixelres(-1)
{
	std::string	sbuf;
	long	i,j,npoly,maxpix=-1, ncap, pixel, ipolyid, icap;
	double weight;
	ifstream	fs(fname.c_str());
	if (!fs) {
		cerr << "Unable to open " << fname << endl;
		exit(1);
	}
	// Read in the number of polygons, which should start the file.
	if (!getline(fs,sbuf)) {cerr<<"Unexpected end-of-file."<<endl;exit(1);}
	istringstream(sbuf) >> npoly;

	double x1,y1,z1,cm1;
	icap = 0;
	while (getline(fs, sbuf)) {
		// Match polygons
		if ((i=sbuf.find("polygon"))!=std::string::npos) {
			// Do all the checks for the different cases here.
			ipolyid = parsepoly(sbuf,ncap,weight,pixel);
			polys.push_back({ipolyid, pixel, ncap, icap, weight});
			icap += ncap;
			if (pixel>maxpix) maxpix=pixel;
			// and read in the caps.
			for (int icap=0; icap<ncap; icap++) {
				if (!getline(fs, sbuf)) throw runtime_error("Unexpected end of file");
				istringstream(sbuf) >> x1 >> y1 >> z1 >> cm1;
				CapClass cap1;
				if (cm1 < 0.0) {
					cap1.x << -1,-x1,-y1,-z1; cap1.cm = cm1;
				} else {
					cap1.x << 1,x1,y1,z1; cap1.cm = cm1;
				}
				caps.push_back(cap1);
			}
			if (caps.size() != icap) throw runtime_error("Egad! Cap and polygon index out of sync");
		}
		// Match pixelization
		else if ((i=sbuf.find("pixelization"))!=std::string::npos &&
				(j=sbuf.find("s",i)         )!=std::string::npos) {
			string ss;
			ss.assign(sbuf,i+string("pixelization").size(),j);
			istringstream(ss) >> pixelres;
			pixeltype = 's';
			cout << "Setting the pixel resolution to " << pixelres << pixeltype << endl;
		}
	}
	fs.close();

	// Check to see that we read in the correct number of polygons
	if (polys.size() != npoly) throw runtime_error("Incorrect number of polygons read in");

	// If we're pixelized make a list of which polygons lie in each pixel.
	if (pixelres>=0) {
		pixelIndex.clear();
		pixelIndex.resize(maxpix+1,-1);
		numPolyInPixel.clear();
		numPolyInPixel.resize(maxpix+1,0);

		// Sort polygon array based on pixelid
		sort(polys.begin(), polys.end(), [](const PolygonClass &p1, const PolygonClass &p2){return p1.pixelid < p2.pixelid;});
		// The cap order is all foobar, make a new copy here.
		CapList tmp = caps;
		icap = 0;
		auto icaptr = caps.begin();
		long ipoly = 0;
		for (PolygonClass &p1 : polys) {
			if (pixelIndex[p1.pixelid] < 0) pixelIndex[p1.pixelid] = ipoly;
			numPolyInPixel[p1.pixelid]++;

			// Copy the caps over to their new location
			icaptr = copy_n(tmp.begin()+p1.capbegin, p1.ncaps, icaptr);
			p1.capbegin = icap;
			icap += p1.ncaps;
			ipoly++;
		}

		// Some runtime checks
		{
			// Check to see number of polygons is correct
			long npolys=0;
			for (long npoly1 : numPolyInPixel) npolys+=npoly1;
			if (npolys != polys.size()) throw runtime_error("Check of numPolyInPixel failed");
		}
		{
			int ipos = 0;
			for (int ndx=0; ndx < pixelIndex.size(); ++ndx) {
				if (pixelIndex[ndx] >= 0) {
					//cout << ndx << " " << pixelIndex[ndx] << " "<< numPolyInPixel[ndx] << endl;
					if (ipos != pixelIndex[ndx]) throw runtime_error("First check of pixelIndex failed");
					for (int jj=0; jj < numPolyInPixel[ndx]; ++jj) {
						if (polys.at(ipos).pixelid != ndx) throw runtime_error("Second check of pixelIndex failed");
						ipos++;
					}
				}
			}
		}

	}


}




