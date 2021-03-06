/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 4)
      cout << "Usage: " << "infile2, which to remove, outfile" << endl;

  ifstream is2(argv[1]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");

  int remove_idx = atoi(argv[2]);

  ofstream of(argv[3]);

  vector<shared_ptr<ftVolume> > volumes;
  
  while (!is2.eof())
    {
      // Read volume from file
      ObjectHeader head;
      is2 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(is2);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      Utils::eatwhite(is2);
    }

  double gap = 0.000001;
  double neighbour = 0.00001;
  double kink = 0.01;
  double bend = 0.05;
  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));


  shared_ptr<ftVolume> remove = model->getBody(remove_idx);
  shared_ptr<ParamVolume> vol2 = remove->getVolume();
  model->removeSolid(remove);
  shared_ptr<ftVolume> vol3 = shared_ptr<ftVolume>(new ftVolume(vol2));
  model->append(vol3);

  int nmb = model->nmbEntities();

  for (int ki=0; ki<nmb; ++ki)
    {
      shared_ptr<ParamVolume> curr = model->getVolume(ki);

      curr->writeStandardHeader(of);
      curr->write(of);
    }
}
