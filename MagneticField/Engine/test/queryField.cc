/** \file
 *  A simple program to print field value.
 *
 *  \author N. Amapane - CERN
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

//#include "DataFormats/GeometryVector/interface/Pi.h"
//#include "DataFormats/GeometryVector/interface/CoordinateSets.h"


#include <iostream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace edm;
using namespace Geom;
using namespace std;

class queryField : public edm::EDAnalyzer {
 public:
  queryField(const edm::ParameterSet& pset) {    
  }

  ~queryField(){}

  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup) {
   ESHandle<MagneticField> magfield;
   setup.get<IdealMagneticFieldRecord>().get(magfield);

   field = magfield.product();

   cout << fixed << setprecision(10);
   cout << "Field Nominal Value: " << field->nominalValue() << endl;

   double x,y,z;

   while (1) {
     
     cout << "Enter X Y Z (cm): ";
    // cout << "Enter R Phi Z (cm): ";
    if (!(cin >> x >>  y >>  z)) exit(0);

      GlobalPoint g(x,y,z);
    // GlobalPoint g(GlobalPoint::Cylindrical(x,y,z));

      GlobalPoint gx(x+0.1,y,z);
      GlobalPoint gy(x,y+0.1,z);
      GlobalPoint gz(x,y,z+0.1);

      GlobalVector b = field->inTesla(g);

      const double bmag = b.mag() ;

      const double bx = field->inTesla(gx).mag();
      const double by = field->inTesla(gy).mag();
      const double bz = field->inTesla(gz).mag();

      const double berr = std::sqrt(std::pow(bmag-bx,2) + std::pow(bmag-by,2) + std::pow(bmag-bz,2));

      cout << "At R=" << g.perp() << " phi=" << g.phi()<< " B=" << b << " |B|=" << b.mag() << " dB=" << berr << endl;// << " BR " << b.transverse() << " BPhi " << b.phi() << endl;
   }
   
  }
   
 private:
  const MagneticField* field;
};


DEFINE_FWK_MODULE(queryField);

