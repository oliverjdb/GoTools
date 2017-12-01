#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "GoParametricTessellableVolume.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"

using namespace Go;
using namespace std;
using namespace TessellateUtils;

int main() {

  // Reading ftVolume
  //string filename = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  string filename = "data/trimmed_cube_graded.g22";
  VolumeModelFileHandler filehandler;
  shared_ptr<ftVolume> testVolume = filehandler.readVolume(filename.c_str());

  // Making a tessellable volume
  GoParametricTessellableVolume ptvolume(*testVolume);

  // Defining target distance between mesh nodes
  const double vdist = 0.04;  
  ostringstream vdist_ss;
  vdist_ss << vdist;

  // Tessellate the volume
  ptvolume.tessellateWithMaterial(vdist);

  // Write the outer shell (triangles) to file
  string plyfile = string("tessellated_shell_")+vdist_ss.str()+string(".ply");
  ofstream os(plyfile);
  ptvolume.writeTessellatedShellPLY(os);
  os.close();

  // Write all points and tetrahedrons to file
  string tetfile = string("tessellated_volume_and_material_")+vdist_ss.str()+string(".tet");
  ofstream os2(tetfile);
  ptvolume.writeTessellatedVolumeAndMaterial(os2);
  os2.close();

  string mshfile = string("tessellated_volume_")+vdist_ss.str()+string(".msh");
  ofstream os3(mshfile);
  ptvolume.writeTessellatedVolumeMSH(os3);
  os3.close();
  
  return 0;
};
  
