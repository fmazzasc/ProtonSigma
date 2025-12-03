
// R__LOAD_LIBRARY(EvtGen)
// R__ADD_INCLUDE_PATH($EVTGEN_ROOT/include)

#include "Pythia8/Pythia.h"
// #include "Pythia8Plugins/ColourReconnectionHooks.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

// #include "Pythia8Plugins/EvtGen.h"
// #include "EvtGen/EvtGen.hh"

// #include "fastjet/PseudoJet.hh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
// #include "fastjet/ClusterSequence.hh"
// #include "fastjet/ClusterSequenceArea.hh"
#include <cstdio> // needed for io
#include <ctime>
#include <iostream> // needed for io
#include <time.h>   /* time */
#include <valarray>
// #include <yaml.h>
//  include <stdio.h>
//  include <glib.h>
// #include <yaml-cpp/yaml.h>

using namespace Pythia8;

int SigmaProton()
{

  // Read JOBID from environment to set unique random seeds per job
  int jobid = 0;
  if (const char *env_p = std::getenv("JOBID"))
  {
    jobid = std::atoi(env_p);
  }

  // Derive a unique seed for Pythia
  int uniqueSeed = 101 + jobid; // 101 is your old seed
  std::cout << "Using jobid = " << jobid << ", seed = " << uniqueSeed << std::endl;

  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Tune:pp = 14");
  // pythia.readString("HardQCD:all = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 13600");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(uniqueSeed));
  pythia.init();
  /*EvtGenDecays *evtgen = 0;
  evtgen = new EvtGenDecays(&pythia, "./DECAY_2010.DEC", "./evt.pdl");
  evtgen->readDecayFile("./DstarDecay.dec");*/

  TTree *fTreeEvent;
  Int_t TreemultSigma;
  Float_t TreeSigmaPx[1000];
  Float_t TreeSigmaPy[1000];
  Float_t TreeSigmaPz[1000];
  Float_t TreeSigmaMass[1000];
  Int_t TreeSigmaMother1[10000];
  Int_t TreeSigmaMother2[10000];
  Int_t TreeSigmaMother1ID[1000];
  Int_t TreeSigmaMother2ID[1000];
  Int_t TreechargeSigma[1000];
  Int_t TreemultProton;
  Float_t TreechargeProton[1000];
  Float_t TreeProtonPx[1000];
  Float_t TreeProtonPy[1000];
  Float_t TreeProtonPz[1000];
  Int_t TreeProtonMother1[10000];
  Int_t TreeProtonMother2[10000];
  Int_t TreeProtonMother1ID[1000];
  Int_t TreeProtonMother2ID[1000];
  fTreeEvent = new TTree("fTreeEvent", "Event");
  fTreeEvent->Branch("TreemultSigma", &TreemultSigma, "TreemultSigma/I");
  fTreeEvent->Branch("TreeSigmaMother1", &TreeSigmaMother1, "TreeSigmaMother1[TreemultSigma]/I");
  fTreeEvent->Branch("TreeSigmaMother2", &TreeSigmaMother2, "TreeSigmaMother2[TreemultSigma]/I");
  fTreeEvent->Branch("TreeSigmaMother1ID", &TreeSigmaMother1ID, "TreeSigmaMother1ID[TreemultSigma]/I");
  fTreeEvent->Branch("TreeSigmaMother2ID", &TreeSigmaMother2ID, "TreeSigmaMother2ID[TreemultSigma]/I");
  fTreeEvent->Branch("TreeSigmaPx", &TreeSigmaPx, "TreeSigmaPx[TreemultSigma]/F");
  fTreeEvent->Branch("TreeSigmaPy", &TreeSigmaPy, "TreeSigmaPy[TreemultSigma]/F");
  fTreeEvent->Branch("TreeSigmaPz", &TreeSigmaPz, "TreeSigmaPz[TreemultSigma]/F");
  fTreeEvent->Branch("TreeSigmaMass", &TreeSigmaMass, "TreeSigmaMass[TreemultSigma]/F");
  fTreeEvent->Branch("TreechargeSigma", &TreechargeSigma, "TreechargeSigma[TreemultSigma]/I");

  fTreeEvent->Branch("TreemultProton", &TreemultProton, "TreemultProton/I");
  fTreeEvent->Branch("TreechargeProton", &TreechargeProton, "TreechargeProton[TreemultProton]/F");
  fTreeEvent->Branch("TreeProtonPx", &TreeProtonPx, "TreeProtonPx[TreemultProton]/F");
  fTreeEvent->Branch("TreeProtonPy", &TreeProtonPy, "TreeProtonPy[TreemultProton]/F");
  fTreeEvent->Branch("TreeProtonPz", &TreeProtonPz, "TreeProtonPz[TreemultProton]/F");
  fTreeEvent->Branch("TreeProtonMother1", &TreeProtonMother1, "TreeProtonMother1[TreemultProton]/I");
  fTreeEvent->Branch("TreeProtonMother2", &TreeProtonMother2, "TreeProtonMother2[TreemultProton]/I");
  fTreeEvent->Branch("TreeProtonMother1ID", &TreeProtonMother1ID, "TreeProtonMother1ID[TreemultProton]/I");
  fTreeEvent->Branch("TreeProtonMother2ID", &TreeProtonMother2ID, "TreeProtonMother2ID[TreemultProton]/I");

  for (int iEvent = 0; iEvent < 3000000; ++iEvent)
  {
    int sigmaNumber = 0;
    int protonNumber = 0;
    if (!pythia.next())
    {
      cout << "check wrong" << "\n";
      continue;
    }
    // evtgen->decay();
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      Particle &part = pythia.event[i];
      Int_t ist = part.status();
      Int_t pid = TMath::Abs(part.id());
      if (std::fabs(part.eta()) > 0.8)
        continue;

      auto mother1 = part.mother1();
      auto mother2 = part.mother2();

      int motherId1 = -999;
      int motherId2 = -999;
      if (mother1 > 0)
      {
        motherId1 = pythia.event[mother1].id();
      }
      if (mother2 >= 0)
      {
        motherId2 = pythia.event[mother2].id();
      }

      if (pid == 3112)
      {
        TreechargeSigma[sigmaNumber] = part.id();
        TreeSigmaPx[sigmaNumber] = part.px();
        TreeSigmaPy[sigmaNumber] = part.py();
        TreeSigmaPz[sigmaNumber] = part.pz();
        TreeSigmaMass[sigmaNumber] = part.m();
        TreeSigmaMother1[sigmaNumber] = part.mother1();
        TreeSigmaMother2[sigmaNumber] = part.mother2();
        TreeSigmaMother1ID[sigmaNumber] = motherId1;
        TreeSigmaMother2ID[sigmaNumber] = motherId2;
        sigmaNumber = sigmaNumber + 1;
      }
      if (pid == 2212)
      {
        TreechargeProton[protonNumber] = part.id();
        TreeProtonPx[protonNumber] = part.px();
        TreeProtonPy[protonNumber] = part.py();
        TreeProtonPz[protonNumber] = part.pz();
        TreeProtonMother1[protonNumber] = part.mother1();
        TreeProtonMother2[protonNumber] = part.mother2();
        TreeProtonMother1ID[protonNumber] = motherId1;
        TreeProtonMother2ID[protonNumber] = motherId2;
        protonNumber = protonNumber + 1;
      }
    }
    if (sigmaNumber * protonNumber > 0)
    {
      TreemultSigma = sigmaNumber;
      TreemultProton = protonNumber;
      fTreeEvent->Fill();
    }
  }

  TFile *fout = new TFile("sigmaprotontree_test.root", "recreate");
  fout->cd();
  fTreeEvent->Write();
  fout->Close();
  pythia.stat();
  return 0;
}
