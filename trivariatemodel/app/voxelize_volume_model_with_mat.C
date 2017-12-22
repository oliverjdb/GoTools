#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include "GoTools/trivariatemodel/lodepng.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
using namespace Go;
using namespace std;
using namespace std::chrono;

int main() {

  // Reading ftVolume
  //string geomfile = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  const char* ifile = "../tessellation/data/nugear_vol_model_with_mat.g22";
  VolumeModelFileHandler fh;
  shared_ptr<VolumeModel> vm = fh.readVolumeModel(ifile);

  int print_dir = 1; // the print direction is y

 
  BoundingBox bb = vm->boundingBox();
  cout << bb << endl;


  // Image (x,y) and layer resolution
  unsigned int resx;
  unsigned int resy;
  unsigned int layerres;
  double stepx;
  double stepy;
  double layerstep; 
  Point vgstep = bb.high()-bb.low();

#if 1 //GEO_RES_SPECIFIED
  stepx = 0.04233; // Stratasys specs
  stepy = 0.08466; // Stratasys specs
  layerstep = 0.027; // Stratasys specs
  layerres = (unsigned int) (bb.high()[print_dir]-bb.low()[print_dir]) / layerstep;
  resx = (unsigned int) (bb.high()[(print_dir+1)%3]-bb.low()[(print_dir+1)%3]) / stepx;
  resy = (unsigned int) (bb.high()[(print_dir+2)%3]-bb.low()[(print_dir+2)%3]) / stepy;
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
  
  cout << "\nStep size: "  << vgstep << "\n" << endl;
  cout << "\nSlice resolution: " << resx << " x " << resy << "\n" << endl;
  cout << "\nNumber of layers: " << layerres << "\n" << endl;

  // Initialize starting point 
  Point vglow(3), low = bb.low();
  vglow[0] = low[0] + vgstep[0]*0.5;
  vglow[1] = low[1] + vgstep[1]*0.5;
  vglow[2] = low[2] + vgstep[2]*0.5;
 
  // Image container
  std::vector<unsigned char> im(resx*resy*4);
  std::vector<int> last_ix(resx*resy,-1);
  std::vector<double> last_uvw(resx*resy*3,0.0);

  double eps = 1.0e-8;

  // Loop through layers
  for (int iz=0;iz!=layerres;++iz) {
    chrono::high_resolution_clock::time_point layerin = high_resolution_clock::now();
    cout << "Start layer" << endl;
    // Set output filename
    stringstream ss;
    ss << iz; //setfill('0') << setw(4) << iz;
    string slayername = string("GearResLayer_")+ss.str()+string(".png");
    const char* layername = slayername.c_str();
    cout << layername << endl;
    unsigned error = lodepng::encode(layername, im, resx, resy);
    // Loop through pixels
    for (int iy=0;iy<resy;++iy) {
      for (int ix=0;ix!=resx;++ix) {
        double u,v,w,d;
        Point vpt;
        if      (print_dir == 2) vpt = vglow+Point(ix*vgstep[0],iy*vgstep[1],iz*vgstep[2]);
        else if (print_dir == 1) vpt = vglow+Point(iy*vgstep[0],iz*vgstep[1],ix*vgstep[2]);
        else if (print_dir == 0) vpt = vglow+Point(iz*vgstep[0],ix*vgstep[1],iy*vgstep[2]);
        //cout << iz << " " << ix << " " << iy << "\n";
        Point clpt;
        int clpidx; 

//cout << "Seed " << last_ix[iy*resx+ix] << " " << last_uvw[3*(iy*resx+ix)] << endl;
//high_resolution_clock::time_point clpin = high_resolution_clock::now();
        vm->closestPoint(vpt,clpidx,u,v,w,clpt,d,eps,last_ix[iy*resx+ix],&(last_uvw[3*(iy*resx+ix)])); 
//high_resolution_clock::time_point clpout = high_resolution_clock::now();
//if (clpidx == last_ix[iy*resx+ix]) cout << "Good seed:" << duration_cast<duration<double>>(clpout-clpin).count() << endl;
//else cout << "Bad seed: " << duration_cast<duration<double>>(clpout-clpin).count() << endl;
        last_ix[iy*resx+ix] = clpidx;
        last_uvw[3*(iy*resx+ix)] = u;
        last_uvw[3*(iy*resx+ix)+1] = v;
        last_uvw[3*(iy*resx+ix)+2] = w; 
        if (fabs(d) < eps && clpidx != -1) {
           shared_ptr<ftVolume> testVolume = vm->getBody(clpidx);  
           vector<double> mpt = testVolume->evaluateMaterialDistribution(u,v,w);
           im[4*resx*iy+4*ix+0] = (unsigned char) (512*mpt[0]);
           im[4*resx*iy+4*ix+1] = 0;
           im[4*resx*iy+4*ix+2] = (unsigned char) (255*mpt[1]);
           im[4*resx*iy+4*ix+3] = 255;
           //if (mpt[0] < 0.2 ) cout << "O";
           //else cout << "X";
        }
        else {
           //cout << "-";
           im[4*resx*iy+4*ix+0] = 0;
           im[4*resx*iy+4*ix+1] = 0;
           im[4*resx*iy+4*ix+2] = 0;
           im[4*resx*iy+4*ix+3] = 255;
        }  
      } 
      //cout << endl;
    }
    cout << im.size() << endl;
    high_resolution_clock::time_point layerout = high_resolution_clock::now();  
    cout << "endlayer" << iz  << " "  << duration_cast<duration<double>>(layerout-layerin).count() << "s" << endl;
  }
 
  return 0;
};
 
