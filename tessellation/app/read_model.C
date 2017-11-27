#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include <iostream>

using namespace std;
using namespace Go;

int main(int varnum, char** vararg) {

  VolumeModelFileHandler fh;

  cout << "Reading ftVolume" << endl;
  shared_ptr<ftVolume> vol = fh.readVolume(vararg[1]);


  cout << "Reading VolumeModel" << endl;
  shared_ptr<VolumeModel> vmodel = fh.readVolumeModel(vararg[1]);
  
  return 0;
};
