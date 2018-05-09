#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include <random>
#include "GoTools/trivariatemodel/lodepng.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
using namespace Go;
using namespace std;
using namespace std::chrono;


enum MaterialMixType{BLENDING, DITHERING}; 


int main( int argc, char* argv[] ) {
  // Reading ftVolume
  //string geomfile = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  if (argc < 3 )
    {
      cerr << "Usage:  " << argv[0] << " input_volume_model_file.g22 output_prefix (slice id and \".png\" will be appended)" << endl;
      return 1;
    }

  const char* ifile = argv[1];
  VolumeModelFileHandler fh;
  shared_ptr<VolumeModel> vm = fh.readVolumeModel(ifile);

  int print_dir = 1; // the print direction is y

  MaterialMixType mix_type = DITHERING;
 
  BoundingBox bb = vm->boundingBox();
  cerr << bb << endl;
 
  // Random number generator for dithering
  std::uniform_real_distribution<double> ran(0.0, 1.0);
  std::default_random_engine dom;
 
  // Image (x,y) and layer resolution
  unsigned int resx;
  unsigned int resy;
  unsigned int layerres;
  double stepx;
  double stepy;
  double layerstep; 
  Point vgstep = bb.high()-bb.low();

#if 1 //GEO_RES_SPECIFIED
  double res_factor = 1; // used to speed up for debugging
  stepx = 0.04233*res_factor; // Stratasys specs
  stepy = 0.08466*res_factor; // Stratasys specs
  layerstep = 0.027*res_factor; // Stratasys specs
  layerres = (unsigned int) ceil((bb.high()[print_dir]-bb.low()[print_dir]) / layerstep) + 2;
  resx = (unsigned int) ceil((bb.high()[(print_dir+1)%3]-bb.low()[(print_dir+1)%3]) / stepx) + 2;
  resy = (unsigned int) ceil((bb.high()[(print_dir+2)%3]-bb.low()[(print_dir+2)%3]) / stepy) + 2;
  vgstep[(print_dir+1)%3] = stepx;
  vgstep[(print_dir+2)%3] = stepy;
#else
  layerres = 12;
  resx = 5;
  resy = 5;
  stepx = (double) (bb.high()[(print_dir+1)%3]-bb.low()[(print_dir+1)%3]) / resx;
  stepy = (double) (bb.high()[(print_dir+2)%3]-bb.low()[(print_dir+2)%3]) / resy;
  layerstep = (bb.high()[print_dir]-bb.low()[print_dir]) / layerres;
  vgstep[(print_dir+1)%3] /= resx;
  vgstep[(print_dir+2)%3] /= resy;
#endif
  vgstep[(print_dir)%3] = layerstep;
 
  cerr << "\nStep size: "  << vgstep << "\n" << endl;
  cerr << "\nSlice resolution: " << resx << " x " << resy << "\n" << endl;
  cerr << "\nNumber of layers: " << layerres << "\n" << endl;

  // Initialize starting point 
  Point vglow(3);
  vglow[0] = bb.low()[0] - vgstep[0]*0.5;
  vglow[1] = bb.low()[1] - vgstep[1]*0.5;
  vglow[2] = bb.low()[2] - vgstep[2]*0.5;
 
  // Image container
  std::vector<unsigned char> im(resx*resy*4);
  std::vector<int> last_ix(resx*resy,-1);
  std::vector<double> last_uvw(resx*resy*3,0.0);

  double eps = 1.0e-8;

  cout << resx << " " << resy << " " << layerres << "\n";

  // Loop through layers
  for (int iz=0;iz!=layerres;++iz) {
    chrono::high_resolution_clock::time_point layerin = high_resolution_clock::now();
    high_resolution_clock::time_point clpin, clpout;
    double clptotal = 0.0;
    cerr << "Start layer " << iz << endl;
    // Set output filename
    stringstream ss;
    ss << iz; //setfill('0') << setw(4) << iz;
    string slayername = string(argv[2])+ss.str()+string(".png");
    const char* layername = slayername.c_str();
    cerr << layername << endl;
    unsigned error = lodepng::encode(layername, im, resx, resy);
    // Loop through pixels
    for (int iy=0;iy<resy;++iy) {
      for (int ix=0;ix!=resx;++ix) {
        double u,v,w,d;
        Point vpt;
        if      (print_dir == 2) vpt = vglow+Point(ix*vgstep[0],iy*vgstep[1],iz*vgstep[2]);
        else if (print_dir == 1) vpt = vglow+Point(iy*vgstep[0],iz*vgstep[1],ix*vgstep[2]);
        else if (print_dir == 0) vpt = vglow+Point(iz*vgstep[0],ix*vgstep[1],iy*vgstep[2]);
        //cerr << iz << " " << ix << " " << iy << "\n";
        Point clpt;
        int clpidx; 

//cerr << "Seed " << last_ix[iy*resx+ix] << " " << last_uvw[3*(iy*resx+ix)] << endl;
        clpin = high_resolution_clock::now();
        vm->closestPoint(vpt,clpidx,u,v,w,clpt,d,eps,last_ix[iy*resx+ix],&(last_uvw[3*(iy*resx+ix)])); 
        clpout = high_resolution_clock::now();
        clptotal += duration_cast<duration<double>>(clpout-clpin).count();

//if (clpidx == last_ix[iy*resx+ix]) cerr << "Good seed:" << duration_cast<duration<double>>(clpout-clpin).count() << endl;
//else cerr << "Bad seed: " << duration_cast<duration<double>>(clpout-clpin).count() << endl;
        last_ix[iy*resx+ix] = clpidx;
        last_uvw[3*(iy*resx+ix)] = u;
        last_uvw[3*(iy*resx+ix)+1] = v;
        last_uvw[3*(iy*resx+ix)+2] = w; 
      
        im[4*resx*iy+4*ix+0] = 0;
        im[4*resx*iy+4*ix+1] = 0;
        im[4*resx*iy+4*ix+2] = 0;
        im[4*resx*iy+4*ix+3] = 255;
        
        if (fabs(d) < eps && clpidx != -1) {
           //cout << "1\n"; 
           shared_ptr<ftVolume> testVolume = vm->getBody(clpidx);  
           vector<double> mpt = testVolume->evaluateMaterialDistribution(u,v,w);
           if (mix_type == DITHERING) {
             double random = ran(dom);
             if (mpt.size() > 0 && random < mpt[0] ) {                  im[4*resx*iy+4*ix+0] = (unsigned char) (255); cout << "1.0\n";}
             else if (mpt.size() > 1 && random < mpt[0]+mpt[1]) {       im[4*resx*iy+4*ix+1] = (unsigned char) (255); cout << "0.5\n";}
             else if (mpt.size() > 2 && random < mpt[0]+mpt[1]+mpt[2]) im[4*resx*iy+4*ix+2] = (unsigned char) (255);
           }
           else if (mix_type == BLENDING) {
             if (mpt.size() > 0) im[4*resx*iy+4*ix+0] = (unsigned char) (255*mpt[0]);
             if (mpt.size() > 1) im[4*resx*iy+4*ix+1] = (unsigned char) (255*mpt[1]);
             if (mpt.size() > 2) im[4*resx*iy+4*ix+2] = (unsigned char) (255*mpt[2]);
             //cerr << (int) im[4*resx*iy+4*ix+0] << " ";
           } 
           else {
             std::cerr << "MaterialMixType not defined\n"; 
             return -1;
           }
           //if (mpt[0] < 0.2 ) cerr << "O";
           //else cerr << "X";
        }
        else {cout << "0\n"; }
      } 
      //cerr << endl;
    }
    cerr << im.size() << endl;

    cerr << "Total time spent on closest point: " << clptotal << endl;
    high_resolution_clock::time_point layerout = high_resolution_clock::now();  
    cerr << "End layer " << iz << " "  << duration_cast<duration<double>>(layerout-layerin).count() << "s" << endl;
  }
 
  return 0;
};
 
