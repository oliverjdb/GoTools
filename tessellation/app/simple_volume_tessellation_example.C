#include <iostream>
#include <fstream>
#include <string>
#include "GoParametricTessellableVolume.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using namespace std;
using namespace TessellateUtils;

int main() {

  // Reading ftVolume
  string filename = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  //string filename = "data/trimmed_cube_tri.g22";
  VolumeModelFileHandler filehandler;
  shared_ptr<ftVolume> testVolume = filehandler.readVolume(filename.c_str());

  // Making a tessellable volume
  GoParametricTessellableVolume ptvolume(*testVolume);

  // Defining target distance between mesh nodes
  const double vdist = 2.0;  

  // Tessellate the volume
  ptvolume.tessellate(vdist);

  // Write the outer shell (triangles) to file
  ofstream os("tessellated_shell.data");
  ptvolume.writeTessellatedShell(os);
  os.close();

  // Write all points and tetrahedrons to file
  ofstream os2("tessellated_volume.data");
  ptvolume.writeTessellatedVolume(os2);
  os2.close();
  
  return 0;
};
  
