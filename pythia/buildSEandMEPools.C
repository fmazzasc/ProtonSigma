

float getkstar(TLorentzVector part1,
               TLorentzVector part2)
{
  TLorentzVector trackSum = part1 + part2;
  TLorentzVector trackRel = part1 - part2;
  const float beta = trackSum.Beta();
  const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
  const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
  const float betaz = beta * std::cos(trackSum.Theta());
  TLorentzVector PartOneCMS(part1);
  TLorentzVector PartTwoCMS(part2);
  const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
  PartOneCMS = boostPRF(PartOneCMS);
  PartTwoCMS = boostPRF(PartTwoCMS);
  const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
  return 0.5 * trackRelK.P();
}

auto setHistogramRangeToContent = [](TH1* h) {
    int first = -1;
    int last  = -1;

    // loop over bins to find first and last non-zero content bin
    int nbins = h->GetNbinsX();
    for (int b = 1; b <= nbins; ++b) {
        if (h->GetBinContent(b) > 0) {
            if (first < 0) first = b;
            last = b;
        }
    }

    if (first < 0 || last < 0) {
        std::cerr << "Warning: histogram " << h->GetName()
                  << " has no bins with content.\n";
        return;
    }

    double xmin = h->GetXaxis()->GetBinLowEdge(first);
    double xmax = h->GetXaxis()->GetBinUpEdge(last);

    h->GetXaxis()->SetRangeUser(xmin, xmax);
};


struct Candidate
{
  TLorentzVector tvec;
  float charge;
  float mother1;
  float mother2;
  float mother1Pdg;
  float mother2Pdg;
};

void buildSEandMEPools(bool reweightPt = true)
{

  // get file with the pT shapes measured in data
  auto filePtShapes = TFile::Open("pt_shapes.root");
  auto hSigmaPtShape = (TH1D *)filePtShapes->Get("hSigmaMinusSignalCounts");
  auto hProtonPtShape = (TH1D *)filePtShapes->Get("hPt");
  hProtonPtShape->Rebin(5); // match binning to sigma

  // set the range of the histograms to the first and last bin with content
  setHistogramRangeToContent(hSigmaPtShape);
  setHistogramRangeToContent(hProtonPtShape);

  hSigmaPtShape->SetDirectory(0);
  hProtonPtShape->SetDirectory(0);
  float sigmaMinPt = hSigmaPtShape->GetXaxis()->GetBinLowEdge(1);
  float sigmaMaxPt = hSigmaPtShape->GetXaxis()->GetBinUpEdge(hSigmaPtShape->GetNbinsX());
  float protonMinPt = hProtonPtShape->GetXaxis()->GetBinLowEdge(1);
  float protonMaxPt = hProtonPtShape->GetXaxis()->GetBinUpEdge(hProtonPtShape->GetNbinsX());
  auto hGenSigmaPtShape = (TH1D *)hSigmaPtShape->Clone("hGenSigmaPtShape");
  hGenSigmaPtShape->Reset();
  auto hGenProtonPtShape = (TH1D *)hProtonPtShape->Clone("hGenProtonPtShape");
  hGenProtonPtShape->Reset();
  auto sigmaPtWeight = (TH1D *)hSigmaPtShape->Clone("sigmaPtWeight");
  auto protonPtWeight = (TH1D *)hProtonPtShape->Clone("protonPtWeight");

  std::vector<std::vector<Candidate>> sigmas;
  std::vector<std::vector<Candidate>> protons;
  // mass, Sigma pt, proton pt, mother flag, kstar, sigma charge, proton charge
  Int_t bins[7] = {700, 250, 250, 2, 300, 3, 3};
  Double_t xmin[7] = {1.0, 0.0, 0.0, 0.0, 0.0, -1.5, -1.5};
  Double_t xmax[7] = {1.4, 10, 10, 2.0, 3.0, 1.5, 1.5};
  double value[7] = {-10.0, -10.0, -10.0, -5.0, -10.0, 0.0, 0.0};
  THnSparseD hsparse_same("hsparse_same", "hsparse_same", 7, bins, xmin, xmax);
  THnSparseD hsparse_mix("hsparse_mix", "hsparse_mix", 7, bins, xmin, xmax);
  TGenPhaseSpace phasespaceevent;
  double masses[2] = {0.13957, 0.93827}; // pion mass, proton mass

  // --- open Sigma–proton tree file ---
  TFile *f = new TFile("trees/lambdaprotontree.root");
  TTree *t = (TTree *)f->Get("fTreeEvent");

  // --- branches (main tree) ---
  TBranch *b_TreemultSigma = nullptr;
  TBranch *b_TreechargeSigma = nullptr;
  TBranch *b_TreeSigmaPx = nullptr;
  TBranch *b_TreeSigmaPy = nullptr;
  TBranch *b_TreeSigmaPz = nullptr;
  TBranch *b_TreeSigmaMass = nullptr;
  TBranch *b_TreeSigmaMother1 = nullptr;
  TBranch *b_TreeSigmaMother2 = nullptr;
  TBranch *b_TreeSigmaMother1ID = nullptr;
  TBranch *b_TreeSigmaMother2ID = nullptr;
  TBranch *b_TreemultProton = nullptr;
  TBranch *b_TreechargeProton = nullptr;
  TBranch *b_TreeProtonPx = nullptr;
  TBranch *b_TreeProtonPy = nullptr;
  TBranch *b_TreeProtonPz = nullptr;
  TBranch *b_TreeProtonMother1 = nullptr;
  TBranch *b_TreeProtonMother2 = nullptr;
  TBranch *b_TreeProtonMother1ID = nullptr;
  TBranch *b_TreeProtonMother2ID = nullptr;
  // leaf types
  Int_t TreemultSigma;
  Int_t TreechargeSigma[1000];
  Float_t TreeSigmaPx[1000];
  Float_t TreeSigmaPy[1000];
  Float_t TreeSigmaPz[1000];
  Float_t TreeSigmaMass[1000];
  Int_t TreeSigmaMother1[10000];
  Int_t TreeSigmaMother2[10000];
  Int_t TreeSigmaMother1ID[1000];
  Int_t TreeSigmaMother2ID[1000];

  Int_t TreemultProton;
  Float_t TreechargeProton[1000];
  Float_t TreeProtonPx[1000];
  Float_t TreeProtonPy[1000];
  Float_t TreeProtonPz[1000];
  Int_t TreeProtonMother1[10000];
  Int_t TreeProtonMother2[10000];
  Int_t TreeProtonMother1ID[1000];
  Int_t TreeProtonMother2ID[1000];

  // set branch addresses (main tree)
  t->SetBranchAddress("TreemultSigma", &TreemultSigma, &b_TreemultSigma);
  t->SetBranchAddress("TreechargeSigma", TreechargeSigma, &b_TreechargeSigma);
  t->SetBranchAddress("TreeSigmaPx", TreeSigmaPx, &b_TreeSigmaPx);
  t->SetBranchAddress("TreeSigmaPy", TreeSigmaPy, &b_TreeSigmaPy);
  t->SetBranchAddress("TreeSigmaPz", TreeSigmaPz, &b_TreeSigmaPz);
  t->SetBranchAddress("TreeSigmaMass", TreeSigmaMass, &b_TreeSigmaMass);
  t->SetBranchAddress("TreeSigmaMother1", TreeSigmaMother1, &b_TreeSigmaMother1);
  t->SetBranchAddress("TreeSigmaMother2", TreeSigmaMother2, &b_TreeSigmaMother2);
  t->SetBranchAddress("TreeSigmaMother1ID", TreeSigmaMother1ID, &b_TreeSigmaMother1ID);
  t->SetBranchAddress("TreeSigmaMother2ID", TreeSigmaMother2ID, &b_TreeSigmaMother2ID);
  t->SetBranchAddress("TreemultProton", &TreemultProton, &b_TreemultProton);
  t->SetBranchAddress("TreechargeProton", TreechargeProton, &b_TreechargeProton);
  t->SetBranchAddress("TreeProtonPx", TreeProtonPx, &b_TreeProtonPx);
  t->SetBranchAddress("TreeProtonPy", TreeProtonPy, &b_TreeProtonPy);
  t->SetBranchAddress("TreeProtonPz", TreeProtonPz, &b_TreeProtonPz);
  t->SetBranchAddress("TreeProtonMother1", TreeProtonMother1, &b_TreeProtonMother1);
  t->SetBranchAddress("TreeProtonMother2", TreeProtonMother2, &b_TreeProtonMother2);
  t->SetBranchAddress("TreeProtonMother1ID", TreeProtonMother1ID, &b_TreeProtonMother1ID);
  t->SetBranchAddress("TreeProtonMother2ID", TreeProtonMother2ID, &b_TreeProtonMother2ID);

  TLorentzVector proton;

  int nevent = t->GetEntries();

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {

    t->GetEntry(jentry);
    int sigmamult = TreemultSigma;
    int protonmult = TreemultProton; // only used if !mixEvents
    Candidate sigmaCand;
    std::vector<Candidate> sigmaVec;
    for (int i = 0; i < sigmamult; ++i)
    {
      sigmaCand.tvec.SetXYZM(TreeSigmaPx[i], TreeSigmaPy[i], TreeSigmaPz[i], TreeSigmaMass[i]);
      sigmaCand.charge = TreechargeSigma[i];
      sigmaCand.mother1 = TreeSigmaMother1[i];
      sigmaCand.mother2 = TreeSigmaMother2[i];
      sigmaCand.mother1Pdg = TreeSigmaMother1ID[i];
      sigmaCand.mother2Pdg = TreeSigmaMother2ID[i];
      if (reweightPt && sigmaCand.tvec.Pt() > sigmaMinPt && sigmaCand.tvec.Pt() < sigmaMaxPt)
      {
        sigmaVec.push_back(sigmaCand);
        hGenSigmaPtShape->Fill(sigmaCand.tvec.Pt());
      }
      else if (!reweightPt)
        sigmaVec.push_back(sigmaCand);
    }
    sigmas.push_back(sigmaVec);

    Candidate protonCand;
    std::vector<Candidate> protonVec;
    for (int j = 0; j < protonmult; ++j)
    {
      protonCand.tvec.SetXYZM(TreeProtonPx[j], TreeProtonPy[j], TreeProtonPz[j], 0.93827);
      protonCand.charge = TreechargeProton[j];
      protonCand.mother1 = TreeProtonMother1[j];
      protonCand.mother2 = TreeProtonMother2[j];
      protonCand.mother1Pdg = TreeProtonMother1ID[j];
      protonCand.mother2Pdg = TreeProtonMother2ID[j];
      if (reweightPt && protonCand.tvec.Pt() > protonMinPt && protonCand.tvec.Pt() < protonMaxPt)
      {
        protonVec.push_back(protonCand);
        hGenProtonPtShape->Fill(protonCand.tvec.Pt());
      }
      else if (!reweightPt)
        protonVec.push_back(protonCand);
    }
    protons.push_back(protonVec);
  }
  std::cout << "Finished loading trees, now processing pairs, SE" << std::endl;
  std::cout << "Sigma sizes: " << sigmas.size() << ", Proton sizes: " << protons.size() << std::endl;

  if (reweightPt)
  {
    // first divide the generated counts for their bin width to get a shape
    for (int i = 1; i <= hGenSigmaPtShape->GetNbinsX(); ++i)
    {
      float binWidth = hGenSigmaPtShape->GetXaxis()->GetBinWidth(i);
      if (hGenSigmaPtShape->GetBinContent(i) > 0)
        hGenSigmaPtShape->SetBinContent(i, hGenSigmaPtShape->GetBinContent(i) / binWidth);
    }
    for (int i = 1; i <= hGenProtonPtShape->GetNbinsX(); ++i)
    {
      float binWidth = hGenProtonPtShape->GetXaxis()->GetBinWidth(i);
      if (hGenProtonPtShape->GetBinContent(i) > 0)
        hGenProtonPtShape->SetBinContent(i, hGenProtonPtShape->GetBinContent(i) / binWidth);
    }
    sigmaPtWeight->Divide(hSigmaPtShape, hGenSigmaPtShape, 1, 1, "B");
    protonPtWeight->Divide(hProtonPtShape, hGenProtonPtShape, 1, 1, "B");

    // normalize weights to have average weight = 1
    double sigmaAvgWeight = sigmaPtWeight->Integral() / sigmaPtWeight->GetNbinsX();
    double protonAvgWeight = protonPtWeight->Integral() / protonPtWeight->GetNbinsX();
    sigmaPtWeight->Scale(1.0 / sigmaAvgWeight);
    protonPtWeight->Scale(1.0 / protonAvgWeight);


  }

  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {
    auto sigmaVec = sigmas[jentry];
    auto protonVec = protons[jentry];
    for (auto &sigma : sigmaVec)
    {
      for (auto &proton : protonVec)
      {

        // create acceptance effect here if needed
        phasespaceevent.SetDecay(sigma.tvec, 2, masses);
        phasespaceevent.Generate();
        auto pidau1 = phasespaceevent.GetDecay(0); // pion
        // discard if mother or daughter out of eta acceptance
        if (std::abs(pidau1->Eta()) > 0.8 || std::abs(sigma.tvec.Eta()) > 0.8 || std::abs(proton.tvec.Eta()) > 0.8)
          continue;

        float relative_momentum = getkstar(sigma.tvec, proton.tvec);
        // mass, Sigma pt, proton pt, mother flag, kstar, sigma charge, proton charge
        value[0] = sigma.tvec.M();
        value[1] = sigma.tvec.Pt();
        value[2] = proton.tvec.Pt();
        value[3] = sigma.mother1 == proton.mother1 ? 0.5 : 1.5;
        value[4] = relative_momentum;
        value[5] = sigma.charge > 0 ? 1.0 : -1.0;
        value[6] = proton.charge > 0 ? 1.0 : -1.0;
        if (std::abs(proton.mother1Pdg) == 3122 || std::abs(proton.mother1Pdg) == 3222)
        {
          continue;
        }

        if (reweightPt)
        {
          double wSigma = sigmaPtWeight->GetBinContent(sigmaPtWeight->FindBin(sigma.tvec.Pt()));
          double wProton = protonPtWeight->GetBinContent(protonPtWeight->FindBin(proton.tvec.Pt()));
          hsparse_same.Fill(value, wProton * wSigma);
        }
        else
        {
          hsparse_same.Fill(value);
        }
      }
    }
  }

  std::cout << "Finished same-event pairs, now processing pairs, ME" << std::endl;

  const int kMaxMixDepth = 4;
  for (Long64_t jentry = 0; jentry < nevent; ++jentry)
  {
    auto &sigmaVec = sigmas[jentry];
    int nMixPartners = std::min<Long64_t>(kMaxMixDepth, nevent - 1);
    for (int idepth = 1; idepth <= nMixPartners; ++idepth)
    {
      Long64_t mixEntry = (jentry + idepth) % nevent;
      auto &protonVec = protons[mixEntry];
      for (auto &sigma : sigmaVec)
      {
        for (auto &proton : protonVec)
        {
          // create acceptance effect here if needed
          phasespaceevent.SetDecay(sigma.tvec, 2, masses);
          phasespaceevent.Generate();
          auto pidau1 = phasespaceevent.GetDecay(0); // pion

          if (std::abs(pidau1->Eta()) > 0.8 ||
              std::abs(sigma.tvec.Eta()) > 0.8 ||
              std::abs(proton.tvec.Eta()) > 0.8)
            continue;

          float relative_momentum = getkstar(sigma.tvec, proton.tvec);
          // mass, Sigma pt, proton pt, mother flag, kstar, sigma charge, proton charge
          value[0] = sigma.tvec.M();
          value[1] = sigma.tvec.Pt();
          value[2] = proton.tvec.Pt();
          value[3] = 1.5; // mixed event: uncorrelated
          value[4] = relative_momentum;
          value[5] = sigma.charge > 0 ? 1.0 : -1.0;
          value[6] = proton.charge > 0 ? 1.0 : -1.0;

          if (std::abs(proton.mother1Pdg) == 3122 || std::abs(proton.mother1Pdg) == 3222)
          {
            continue;
          }

          if (reweightPt)
          {
            double wSigma = sigmaPtWeight->GetBinContent(sigmaPtWeight->FindBin(sigma.tvec.Pt()));
            double wProton = protonPtWeight->GetBinContent(protonPtWeight->FindBin(proton.tvec.Pt()));
            hsparse_mix.Fill(value, wProton * wSigma);
          }
          else
          {
            hsparse_mix.Fill(value);
          }
        }
      }
    }
  }

  // dump to output file
  TFile *fout = new TFile("lambdaproton_thnsparse_rew.root", "recreate");
  fout->cd();
  hsparse_same.Write();
  hsparse_mix.Write();
  if (reweightPt)
  {

    hGenSigmaPtShape->Write();
    hGenProtonPtShape->Write();
    hProtonPtShape->Write();
    hSigmaPtShape->Write();
    sigmaPtWeight->Write();
    protonPtWeight->Write();
  }
  fout->Close();
}