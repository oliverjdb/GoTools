#include <iostream>
#include <fstream>
#include <string>
#include "GoParametricTesselableVolume.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using namespace std;
using namespace TesselateUtils;

int main() {

  // Reading ftVolume
  string filename = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  VolumeModelFileHandler filehandler;
  shared_ptr<ftVolume> testVolume = filehandler.readVolume(filename.c_str());

  // Making a tesselable volume
  GoParametricTesselableVolume ptvolume(*testVolume);

  // Defining target distance between mesh nodes
  const double vdist = 0.5;  

  // Tesselate the volume
  ptvolume.tesselate(vdist);

  // Write the outer shell (triangles) to file
  ofstream os("tesselated_shell.data");
  ptvolume.writeTesselatedShell(os);
  os.close();

  // Write all points and tetrahedrons to file
  ofstream os2("tesselated_volume.data");
  ptvolume.writeTesselatedVolume(os2);
  os2.close();
  
  return 0;
};
  
