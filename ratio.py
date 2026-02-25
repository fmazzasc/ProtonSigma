import ROOT
import argparse

parser = argparse.ArgumentParser(description='Calculate the ratio of k* distributions from MC and ME data.')
parser.add_argument('--mc', action='store_true', help='Analyze MC data')
parser.add_argument('--bantib', action='store_true', help='Analyze b-anti-b pairs')
parser.add_argument('--output_dir', type=str, default='results_recal', help='Directory to save the output ROOT file')
args = parser.parse_args()

is_mc = args.mc
is_bantib = args.bantib
output_dir = args.output_dir
suffix_bantib = is_bantib and "_bantib" or ""
suffix_mc = "_MC" if is_mc else ""


data_file = ROOT.TFile(f"{output_dir}/kstar_sigma_proton_analysis{suffix_mc}{suffix_bantib}.root")
kstar_histo = data_file.Get("hKstarData")
kstar_histo.SetDirectory(0)
## draw with pe
kstar_histo.SetMarkerStyle(20)
# kstar_histo.Rebin(2)  
em_file = ROOT.TFile(f"{output_dir}/kstar_sigma_proton_analysis{suffix_mc}_ME{suffix_bantib}.root")
kstar_histo_em = em_file.Get("hKstarData")
kstar_histo_em.SetDirectory(0)
kstar_histo_em.SetMarkerStyle(20)
# kstar_histo_em.Rebin(2)
ratio_histo = kstar_histo.Clone("hRatio")
##scale ME to the same integral in 0.4 - 0.6 GeV/c
bin_low = ratio_histo.GetXaxis().FindBin(0.5)
bin_high = ratio_histo.GetXaxis().FindBin(0.8)
integral_data = ratio_histo.Integral(bin_low, bin_high)
integral_em = kstar_histo_em.Integral(bin_low, bin_high)
scaling_factor = integral_data / integral_em
kstar_histo_em.Scale(scaling_factor)
ratio_histo.Divide(kstar_histo_em)

kstar_histo_sideband = data_file.Get("hKStarSideband")
kstar_histo_sideband.SetDirectory(0)
kstar_histo_sideband.SetMarkerStyle(20)
# kstar_histo_sideband.Rebin(2)
kstar_histo_sideband_em = em_file.Get("hKStarSideband")
kstar_histo_sideband_em.SetDirectory(0)
kstar_histo_sideband_em.SetMarkerStyle(20)
# kstar_histo_sideband_em.Rebin(2)
ratio_histo_sideband = kstar_histo_sideband.Clone("hRatioSideband")
integral_data_sb = ratio_histo_sideband.Integral(bin_low, bin_high)
integral_em_sb = kstar_histo_sideband_em.Integral(bin_low, bin_high)
scaling_factor_sb = integral_data_sb / integral_em_sb
kstar_histo_sideband_em.Scale(scaling_factor_sb)
ratio_histo_sideband.Divide(kstar_histo_sideband_em)



if is_mc:
    kstar_histo_signal = data_file.Get("hKstarSignal")
    kstar_histo_signal.SetDirectory(0)
    kstar_histo_bkg = data_file.Get("hKstarBkg")
    kstar_histo_bkg.SetDirectory(0)
    # kstar_histo_signal.Rebin(2)
    # kstar_histo_bkg.Rebin(2)

    kstar_histo_signal_em = em_file.Get("hKstarSignal")
    kstar_histo_signal_em.SetDirectory(0)
    kstar_histo_bkg_em = em_file.Get("hKstarBkg")
    kstar_histo_bkg_em.SetDirectory(0)
    # kstar_histo_signal_em.Rebin(2)
    # kstar_histo_bkg_em.Rebin(2)

    ratio_histo_signal = kstar_histo_signal.Clone("hRatioSignal")
    ratio_histo_bkg = kstar_histo_bkg.Clone("hRatioBkg")
    integral_data_signal = kstar_histo_signal.Integral(bin_low, bin_high)
    integral_em_signal = kstar_histo_signal_em.Integral(bin_low, bin_high)
    scaling_factor_signal = integral_data_signal / integral_em_signal
    kstar_histo_signal_em.Scale(scaling_factor_signal)
    ratio_histo_signal.Divide(kstar_histo_signal_em)
    integral_data_bkg = kstar_histo_bkg.Integral(bin_low, bin_high)
    integral_em_bkg = kstar_histo_bkg_em.Integral(bin_low, bin_high)
    scaling_factor_bkg = integral_data_bkg / integral_em_bkg
    kstar_histo_bkg_em.Scale(scaling_factor_bkg)
    ratio_histo_bkg.Divide(kstar_histo_bkg_em)




outfile = ROOT.TFile(f"{output_dir}/kstar_ratio_sigma_proton{suffix_mc}{suffix_bantib}.root", "RECREATE")
kstar_histo_em.Write()
ratio_histo.Write()
ratio_histo_sideband.Write()
if is_mc:
    ratio_histo_signal.Write()
    ratio_histo_bkg.Write()
outfile.Close()
