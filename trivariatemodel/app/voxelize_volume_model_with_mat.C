#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "GoTools/trivariatemodel/lodepng.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using namespace std;


int main() {

  // Reading ftVolume
  //string geomfile = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  const char* ifile = "../tessellation/data/nugear_vol_model_with_mat.g22";
  VolumeModelFileHandler fh;
  shared_ptr<VolumeModel> vm = fh.readVolumeModel(ifile);

  int print_dir = 0; // the print direction is y

  //if (print_dir ) 
 
  BoundingBox bb = vm->boundingBox();
  cout << bb << endl;

  

  // Image (x,y) and layer resolution
  unsigned int resx=256,resy=256,layerres=256;

  
  // Bottom left pixel
  Point vglow(3), vgstep(3), low = bb.low();
  vgstep = bb.high()-bb.low(); 
  if      (print_dir == 2) {
    vgstep[0] /= resx;
    vgstep[1] /= resy;
    vgstep[2] /= layerres; 
 }
  else if (print_dir == 1) {
    vgstep[2] /= resx;
    vgstep[0] /= resy;
    vgstep[1] /= layerres;
  }
  else if (print_dir == 0) {
    vgstep[1] /= resx;
    vgstep[2] /= resy;
    vgstep[0] /= layerres;
  }
  else return -1;
 
  // Initialize starting point 
  vglow[0] = low[0] + vgstep[0]*0.5;
  vglow[1] = low[1] + vgstep[1]*0.5;
  vglow[2] = low[2] + vgstep[2]*0.5;
 
  // Image container
  std::vector<unsigned char> im(resx*resy*4);

  double eps = 1.0e-8;

  // Loop through layers
  for (int iz=0;iz!=layerres;++iz) {
    cout << "Start layer" << endl;
    // Set output filename
    stringstream ss;
    ss << setfill('0') << setw(4) << iz;
    const char* layername = (string("layer")+ss.str()+string(".png")).c_str();
    cout << layername << endl;
    unsigned error = lodepng::encode(layername, im, resx, resy);
    // Loop through pixels
    for (int iy=0;iy!=resy;++iy) {
      for (int ix=0;ix!=resx;++ix) {
        double u,v,w,d;
        Point vpt;
        if      (print_dir == 2) vpt = vglow+Point(ix*vgstep[0],iy*vgstep[1],iz*vgstep[2]);
        else if (print_dir == 1) vpt = vglow+Point(iz*vgstep[0],ix*vgstep[1],iy*vgstep[2]);
        else if (print_dir == 0) vpt = vglow+Point(iy*vgstep[0],iz*vgstep[1],ix*vgstep[2]);
        //cout << vpt << endl;
        Point clpt;
        int clpidx;
        vm->closestPoint(vpt,clpidx,u,v,w,clpt,d,eps); 
        if (fabs(d) < 0.0000001 && clpidx != -1) {
           shared_ptr<ftVolume> testVolume = vm->getBody(clpidx);  
           vector<double> mpt = testVolume->evaluateMaterialDistribution(u,v,w);
           im[4*resx*iy+4*ix+0] = (unsigned char) (512*mpt[0]);
           im[4*resx*iy+4*ix+1] = 0;
           im[4*resx*iy+4*ix+2] = (unsigned char) (255*mpt[1]);
           im[4*resx*iy+4*ix+3] = 255;
           if (mpt[0] < 0.2 ) cout << "O";
           else cout << "X";
        }
        else {
           cout << "-";
           im[4*resx*iy+4*ix+0] = 0;
           im[4*resx*iy+4*ix+1] = 0;
           im[4*resx*iy+4*ix+2] = 0;
           im[4*resx*iy+4*ix+3] = 255;
        }  
      } 
      cout << endl;
    }
    cout << im.size() << endl;
    cout << "endlayer" << iz  << " "  << error << endl;
  }
 
  return 0;
};
 
