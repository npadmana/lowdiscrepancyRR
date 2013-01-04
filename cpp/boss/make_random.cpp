using namespace std;

#include	<cmath>
#include	<iostream>
#include	<fstream>
#include	<sstream>
#include	<iomanip>
#include	<vector>
#include	<list>
#include	<string>

/*

Mask class for manipulating the survey mask.

Author:	Martin White	(UCB)
Written:	10-Jun-2009
Modified:	 6-Aug-2009	(Put MaskClass in Mangle namespace)

*/

namespace Mangle {

class CapClass {
private:
  double	x,y,z,cm;
public:
  CapClass(double x1, double y1, double z1, double cm1) {
    x = x1; y = y1; z = z1; cm = cm1;
  }
  bool incap(double x0, double y0, double z0) {
    double	cdot;
    bool	tmp;
    cdot = 1.0 - x*x0 - y*y0 - z*z0;
    if (cm < 0.0)
      tmp = cdot > (-cm);
    else
      tmp = cdot < cm;
    return(tmp);
  }
  bool incap(double theta, double phi) {
    double	cdot,x0,y0,z0;
    bool	tmp;
    x0 = sin(theta)*cos(phi);
    y0 = sin(theta)*sin(phi);
    z0 = cos(theta); 
    cdot= 1.0 - x*x0 - y*y0 - z*z0;
    if (cm < 0.0)
      tmp = cdot > (-cm);
    else
      tmp = cdot < cm;
    return(tmp);
  }
};	// CapClass



class PolygonClass {
private:
  std::list<CapClass>	caps;
  long			pixelid;
  double		weight;
public:
  PolygonClass() {
    caps.clear();
    weight = 0;
    pixelid    = 0;
  }
  long getpixelid() {
    return(pixelid);
  }
  void setpixelid(long pid1) {
    pixelid = pid1;
  }
  double getwt() {
    return(weight);
  }
  void setwt(double wt) {
    weight = wt;
  }
  void addcap(CapClass cap1) {
    caps.push_back(cap1);
  }
  bool inpoly(double x, double y, double z) {	// Not used currently.
    bool tmp=true;
    for (std::list<CapClass>::iterator i=caps.begin();
         i!=caps.end() && tmp; i++) {
      tmp = tmp && ( i->incap(x,y,z) );
    }
    return(tmp);
  }
  bool inpoly(double theta, double phi) {
    bool tmp=true;
    for (std::list<CapClass>::iterator i=caps.begin();
         i!=caps.end() && tmp; i++) {
      tmp = tmp && ( i->incap(theta,phi) );
    }
    return(tmp);
  }
};	// PolygonClass


class SDSSpixClass {	// Does conversions for us.  Based on SDSSpix.pro
private:
  long		nx0,ny0,resolution;
  double	deg2Rad,node,etapole;
  double	etaOffSet,surveyCenterRA,surveyCenterDEC;
  std::vector<double> tp2survey(double theta, double phi) {
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
public:
  SDSSpixClass(long pixres) {
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
  long pixelnum(double theta, double phi) {
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
};	// SDSSpixClass




class MaskClass {
// A class to handle masks.  When invoked it reads an ascii mask file.
// The main method returns a completeness, used for making random catalogs.
private:
  std::vector<PolygonClass>	polygons;
  std::vector< std::list<int> >	pixels;
  int				pixelres;
  char				pixeltype;
  long parsepoly(std::string sbuf, long &ncap, double &weight, long &pixel) {
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
  long pixelnum(double theta, double phi) {
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

public:
  MaskClass(std::string fname) {	// Load the mask from fname.
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
      polygons[ipoly].setpixelid(pixel);
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
        pixels[polygons[ipoly].getpixelid()].push_back(ipoly);
    }
  }
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

}	// namespace Mangle





int	main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "Usage: make_random <polygon-file>" << endl;
    exit(1);
  }
  std::string fname(argv[1]);
  Mangle::MaskClass mm(fname);
  cout << "# Read " << mm.npolygons() << " polygons from " << fname << endl;
  cout << setiosflags(ios::left)  << "#" << setiosflags(ios::right)
       << setw(9) << "Phi" << setw(10) << "Cos(th)" << endl;
  int i=0;
  do {
    double cth = 1.5*drand48()-0.5;
    double phi = 2*M_PI*drand48();
    if (mm.completeness(acos(cth),phi)>0.90) {
      cout << fixed << setprecision(5)
           << setw(10) << phi << setw(10) << cth << endl;
      i++;
    }
    if (i%100==1) cout.flush();
  } while (i<5000);
  return(0);
}
