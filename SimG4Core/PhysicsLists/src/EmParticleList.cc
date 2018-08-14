#include "SimG4Core/PhysicsLists/interface/EmParticleList.h"

EmParticleList::EmParticleList() 
{
  pNames = 
    { 
        "gamma",            "e-",           "e+",           "mu+",        "mu-",
          "pi+",           "pi-",        "kaon+",         "kaon-",     "proton",
  "anti_proton",         "alpha",          "He3",    "GenericIon",         "B+",
           "B-",            "D+",           "D-",           "Ds+",        "Ds-",
     "anti_He3",    "anti_alpha","anti_deuteron","anti_lambda_c+","anti_omega-",
"anti_sigma_c+","anti_sigma_c++",  "anti_sigma+",   "anti_sigma-","anti_triton",
   "anti_xi_c+",      "anti_xi-",     "deuteron",     "lambda_c+",     "omega-",
     "sigma_c+",     "sigma_c++",       "sigma+",        "sigma-",       "tau+",
         "tau-",        "triton",        "xi_c+",           "xi-"
    };
}

EmParticleList::~EmParticleList() 
{}

const std::vector<G4String>& EmParticleList::PartNames() const
{
  return pNames;
}

