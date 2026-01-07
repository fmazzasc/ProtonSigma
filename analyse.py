import ROOT
import argparse

ROOT.gInterpreter.Declare("""
float calcKstar(float pxSig, float pySig, float pzSig, float pxPr, float pyPr, float pzPr, float mSig, float mPr) {
    TLorentzVector sigVec;  
    TLorentzVector prVec;
    sigVec.SetXYZM(pxSig, pySig, pzSig, mSig);
    prVec.SetXYZM(pxPr, pyPr, pzPr, mPr);
    TLorentzVector trackSum = sigVec + prVec;
    float beta = trackSum.Beta();
    float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    float betaz = beta * cos(trackSum.Theta());
    TLorentzVector PartOneCMS = sigVec;
    TLorentzVector PartTwoCMS = prVec;
    ROOT::Math::Boost boostPRF(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
}
""")

ROOT.gInterpreter.Declare("""
float calcMass(float px1, float py1, float pz1, float m1, float px2, float py2, float pz2, float m2) {
    TLorentzVector vec1;
    TLorentzVector vec2;
    px1 -= px2;
    py1 -= py2;
    pz1 -= pz2;
    vec1.SetPxPyPzE(px1, py1, pz1, sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1));
    vec2.SetPxPyPzE(px2, py2, pz2, sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2));
    TLorentzVector vecSum = vec1 + vec2;
    return vecSum.M();
}
""")

ROOT.gInterpreter.Declare("""
float calcKinkAngle(float pxMother, float pyMother, float pzMother, float pxDaughter, float pyDaughter, float pzDaughter) {
    TLorentzVector motherVec;
    TLorentzVector daughterVec;
    motherVec.SetPxPyPzE(pxMother, pyMother, pzMother, sqrt(pxMother*pxMother + pyMother*pyMother + pzMother*pzMother + 1.1965*1.1965));
    daughterVec.SetPxPyPzE(pxDaughter, pyDaughter, pzDaughter, sqrt(pxDaughter*pxDaughter + pyDaughter*pyDaughter + pzDaughter*pzDaughter + 0.13957*0.13957));
    return motherVec.Angle(daughterVec.Vect());
}
""")

## enable parallel processing
ROOT.ROOT.EnableImplicitMT()

## add argparse for mc and me options
parser = argparse.ArgumentParser(description='Analyze Sigma-Proton pairs to calculate k* and other variables.')
parser.add_argument('--mc', action='store_true', help='Analyze MC data')
parser.add_argument('--me', action='store_true', help='Analyze ME data')
parser.add_argument('--bantib', action='store_true', help='Analyze b-anti-b pairs')

args = parser.parse_args()
is_mc = args.mc
is_me = args.me
is_bantib = args.bantib
bantib_cut = is_bantib and "fChargeSigma * fChargePr > 0" or "fChargeSigma * fChargePr < 0"
suffix_bantib = is_bantib and "_bantib" or ""
print('----------------------------------')
print("Analyzing with options - MC:", is_mc, " ME:", is_me, " b-anti-b:", is_bantib)

suffix_mc = "_MC" if is_mc else ""
suffix_me = "_ME" if is_me else ""
mc_dir = "mc"
data_dir = "23_thin"
dir_to_use = mc_dir if is_mc else data_dir

input_data = f"{dir_to_use}/AO2D{suffix_mc}{suffix_me}.root"

tree_name = is_mc and "O2sigmaprotonmcca" or "O2sigmaprotoncand"
file_data_list = input_data if isinstance(input_data, list) else [input_data]
chainData = ROOT.TChain("O2sigmaprotoncand")
for fileName in file_data_list:
  fileData = ROOT.TFile(fileName)
  for key in fileData.GetListOfKeys() :
    keyName = key.GetName()
    if 'DF_' in keyName :
        chainData.Add(f'{fileName}/{keyName}/{tree_name}')
dataDf = ROOT.RDataFrame(chainData)
print(f"Data entries: {dataDf.Count().GetValue()}")
## calc kstar and add as new columns


dataDf = dataDf.Define("fKstar", "calcKstar(fPxMoth, fPyMoth, fPzMoth, fPxPr, fPyPr, fPzPr, 1.1965, 0.93827)")  # masses in GeV/c2
dataDf = dataDf.Define("fMassSigma", "calcMass(fPxMoth, fPyMoth, fPzMoth, 0.938272, fPxDaug, fPyDaug, fPzDaug, 0.13957)")  # masses in GeV/c2
dataDf = dataDf.Define("fKinkAngle", "calcKinkAngle(fPxMoth, fPyMoth, fPzMoth, fPxDaug, fPyDaug, fPzDaug)")  # masses in GeV/c2
dataDf = dataDf.Define("fPtProton", "sqrt(fPxPr * fPxPr + fPyPr * fPyPr)")
dataDf = dataDf.Define("fSigmaPt", "sqrt(fPxMoth * fPxMoth + fPyMoth * fPyMoth)")



dataDf = dataDf.Filter(bantib_cut)
dataSideBand = dataDf.Filter("fMassSigma > 1.25")  # side-band upper
dataDf = dataDf.Filter("fSigmaPt > 1.2 && fSigmaPt < 3")  # pt cut
hMassVsKStarBefCuts = dataDf.Histo2D(("hMassVsKStarBefCuts", ";Mass #Sigma (GeV/#it{c}^{2});k* (GeV/#it{c})", 100, 1.1, 1.3, 100, 0., 3), "fMassSigma", "fKstar")
dataDf = dataDf.Filter("fMassSigma > 1.19 && fMassSigma < 1.21")  # mass window cut

## create histogram of kstar
hKstarData = dataDf.Histo1D(("hKstarData", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar")
hMassSigmaData = dataDf.Histo1D(("hMassSigmaData", ";Mass #Sigma (GeV/#it{c}^{2});Counts", 100, 1.1, 1.3), "fMassSigma")
hKinkAngleData = dataDf.Histo1D(("hKinkAngleData", ";Kink Angle (rad);Counts", 100, 0., 1.5), "fKinkAngle")
hPtProtonData = dataDf.Histo1D(("hPtProtonData", ";#it{p}_{T} proton (GeV/#it{c});Counts", 100, 0., 5), "fPtProton")
hPtSigmaData = dataDf.Histo1D(("hPtSigmaData", ";#it{p}_{T} #Sigma (GeV/#it{c});Counts", 100, 0., 5), "fSigmaPt")
hKStarVsMassSigma = dataDf.Histo2D(("hKStarVsMassSigma", ";Mass #Sigma (GeV/#it{c}^{2});k* (GeV/#it{c})", 100, 1.1, 1.3, 100, 0., 3), "fMassSigma", "fKstar")
hKStarVsPtSigma = dataDf.Histo2D(("hKStarVsPtSigma", ";#it{p}_{T} #Sigma (GeV/#it{c});k* (GeV/#it{c})", 100, 0., 5, 100, 0., 3), "fSigmaPt", "fKstar")
hKStarSideband = dataSideBand.Histo1D(("hKStarSideband", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar")

if is_mc:
    dataDfSigma = dataDf.Filter("abs(fSigmaPDG)==3112")  # signal sample
    dataDfSigmaSideband = dataSideBand.Filter("abs(fSigmaPDG)==3112")  # signal sample in sideband
    dataBkg = dataDf.Filter("abs(fSigmaPDG)!=3112")  # background sample
    dataSigmaPlus = dataDf.Filter("abs(fSigmaPDG)==3222")
    ## create histogram of kstar for background

    hDecRadVsKstarSigmaMinus = dataDfSigma.Histo2D(("hDecRadVsKstarSigmaMinus", ";Decay Radius (cm);k* (GeV/#it{c})", 100, 0., 300., 100, 0., 3), "fSigmaDecRad", "fKstar")
    hDecRadVsKstarSigmaPlus = dataSigmaPlus.Histo2D(("hDecRadVsKstarSigmaPlus", ";Decay Radius (cm);k* (GeV/#it{c})", 100, 0., 300., 100, 0., 3), "fSigmaDecRad", "fKstar")

    hKstarSignal = dataDfSigma.Histo1D(("hKstarSignal", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar")
    hKstarSignalSideband = dataDfSigmaSideband.Histo1D(("hKstarSignalSideband", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar")
    hMassSigmaSignal = dataDfSigma.Histo1D(("hMassSigmaSignal", ";Mass #Sigma (GeV/#it{c}^{2});Counts", 100, 1.1, 1.3), "fMassSigma")
    hKinkAngleSignal = dataDfSigma.Histo1D(("hKinkAngleSignal", ";Kink Angle (rad);Counts", 100, 0., 1.5), "fKinkAngle")
    hPtSignal = dataDfSigma.Histo1D(("hPtSigmaSignal", ";#it{p}_{T} #Sigma (GeV/#it{c});Counts", 100, 0., 5), "fSigmaPt")

    hKstarBkg = dataBkg.Histo1D(("hKstarBkg", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar") 
    hMassSigmaBkg = dataBkg.Histo1D(("hMassSigmaBkg", ";Mass #Sigma (GeV/#it{c}^{2});Counts", 100, 1.1, 1.3), "fMassSigma")
    hKinkAngleBkg = dataBkg.Histo1D(("hKinkAngleBkg", ";Kink Angle (rad);Counts", 100, 0., 1.5), "fKinkAngle")

    hPurityKstar = hKstarSignal.Clone("hPurityKstar")
    hPurityKstar.Divide(hKstarData.GetPtr())

    hPurityKstarSideband = hKstarSignalSideband.Clone("hPurityKstarSideband")
    hPurityKstarSideband.Divide(hKStarSideband.GetPtr())

    hPurityPtSigma = hPtSignal.Clone("hPurityPtSigma")
    hPurityPtSigma.Divide(hPtSigmaData.GetPtr())

    hSigmaPlusKstar = dataSigmaPlus.Histo1D(("hSigmaPlusKstar", ";k* (GeV/#it{c});Counts", 100, 0., 3), "fKstar")   
    hSigmaPlusToMinusRatio = hSigmaPlusKstar.Clone("hSigmaPlusToMinusRatio")
    hSigmaPlusToMinusRatio.Divide(hKstarSignal.GetPtr())
    hSigmaPlusDauPdg = dataSigmaPlus.Histo1D(("hSigmaPlusDauPdg", ";Daughter PDG;Counts", 2000, -1000, 1000), "fDaughterPDG")
    dataBkgLowKstar = dataBkg.Filter("fKstar < 0.2")
    hSigmaPlusKStarPt = dataSigmaPlus.Histo2D(("hSigmaPlusKStarPt", ";#it{p}_{T} #Sigma (GeV/#it{c});k* (GeV/#it{c})", 100, 0., 5, 100, 0., 3), "fSigmaPt", "fKstar")
    hSigmaMinusKStarPt = dataDfSigma.Histo2D(("hSigmaMinusKStarPt", ";#it{p}_{T} #Sigma (GeV/#it{c});k* (GeV/#it{c})", 100, 0., 5, 100, 0., 3), "fSigmaPt", "fKstar")

    hPDGMom = dataBkgLowKstar.Histo1D(("hPDGMom", ";Mother PDG;Counts", 2000, -10000, 10000), "fSigmaPDG")
    hPDGDau = dataBkgLowKstar.Histo1D(("hPDGDau", ";Daughter PDG;Counts", 2000, -10000, 10000), "fDaughterPDG")
    # hPDGMomVsDau = dataBkgLowKstar.Histo2D(("hPDGMomVsDau", ";Mother PDG;Daughter PDG", 20000, -10000, 10000, 2000, -10000, 10000), "fSigmaPDG", "fDaughterPDG")
    

outfile = ROOT.TFile(f"results/kstar_sigma_proton_analysis{suffix_mc}{suffix_me}{suffix_bantib}.root", "recreate")
hMassVsKStarBefCuts.Write()
hKstarData.Write()
hMassSigmaData.Write()
hKinkAngleData.Write()
hPtProtonData.Write()
hPtSigmaData.Write()
hKStarVsMassSigma.Write()
hKStarVsPtSigma.Write()
hKStarSideband.Write()


if is_mc:
    hKstarSignal.Write()
    hMassSigmaSignal.Write()
    hKinkAngleSignal.Write()
    hDecRadVsKstarSigmaMinus.Write()
    hDecRadVsKstarSigmaPlus.Write()
    hKstarBkg.Write()
    hMassSigmaBkg.Write()
    hKinkAngleBkg.Write()
    hPurityKstar.Write()
    hPurityPtSigma.Write()
    hPurityKstarSideband.Write()
    hSigmaPlusKstar.Write()
    hSigmaPlusToMinusRatio.Write()
    hSigmaPlusDauPdg.Write()
    hPDGMom.Write()
    hPDGDau.Write()
    # hPDGMomVsDau.Write()

    hSigmaPlusKStarPt.Write()
    hSigmaMinusKStarPt.Write()



outfile.Close()