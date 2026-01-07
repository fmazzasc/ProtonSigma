import ROOT

n_mix = 4

file_input = ROOT.TFile.Open("sigmaproton_thnsparse_rew.root")
thnsparse_same = file_input.Get("hsparse_same")
thnsparse_mix = file_input.Get("hsparse_mix")

out_file = ROOT.TFile.Open("kstar_ratio_sigmaproton_rew.root", "RECREATE")

# values: mass, Sigma pt, proton pt, mother flag, kstar, sigma charge, proton charge
baryon_types = {'baryon_baryon':[1,1], 'barion_antibaryon':[1, -1],'antibaryon_antibaryon':[-1,-1]}
h_sigma_pt = thnsparse_same.Projection(1)
h_proton_pt = thnsparse_same.Projection(2)

for baryon_type, baryon_charges in baryon_types.items():
    ## reset the ranges
    thnsparse_same.GetAxis(3).SetRange(1, thnsparse_same.GetAxis(3).GetNbins())
    thnsparse_mix.GetAxis(3).SetRange(1, thnsparse_mix.GetAxis(3).GetNbins())
    thnsparse_same.GetAxis(5).SetRange(1, thnsparse_same.GetAxis(5).GetNbins())
    thnsparse_mix.GetAxis(5).SetRange(1, thnsparse_mix.GetAxis(5).GetNbins())
    thnsparse_same.GetAxis(6).SetRange(1, thnsparse_same.GetAxis(6).GetNbins())
    thnsparse_mix.GetAxis(6).SetRange(1, thnsparse_mix.GetAxis(6).GetNbins())
    thnsparse_same.GetAxis(5).SetRange(thnsparse_same.GetAxis(5).FindBin(baryon_charges[0]), thnsparse_same.GetAxis(5).FindBin(baryon_charges[0]))
    thnsparse_same.GetAxis(6).SetRange(thnsparse_same.GetAxis(6).FindBin(baryon_charges[1]), thnsparse_same.GetAxis(6).FindBin(baryon_charges[1]))
    thnsparse_mix.GetAxis(5).SetRange(thnsparse_mix.GetAxis(5).FindBin(baryon_charges[0]), thnsparse_mix.GetAxis(5).FindBin(baryon_charges[0]))
    thnsparse_mix.GetAxis(6).SetRange(thnsparse_mix.GetAxis(6).FindBin(baryon_charges[1]), thnsparse_mix.GetAxis(6).FindBin(baryon_charges[1]))
    out_file.mkdir(baryon_type)
    out_file.cd(baryon_type)

    k_star_ratio_inclusive = thnsparse_same.Projection(4).Clone(f"k_star_ratio_inclusive_{baryon_type}")
    k_star_ratio_inclusive.Sumw2()
    k_star_ratio_inclusive.Rebin(10)
    k_star_ratio_inclusive.GetXaxis().SetTitle("#it{k*} (GeV/#it{c})")
    k_star_ratio_inclusive.GetYaxis().SetTitle("C(#it{k*})")

    thnsparse_same.GetAxis(3).SetRange(thnsparse_same.GetAxis(3).FindBin(1.5), thnsparse_same.GetAxis(3).FindBin(1.5))
    k_star_ratio_diff_mother = thnsparse_same.Projection(4).Clone(f"k_star_ratio_diff_mother_{baryon_type}")
    k_star_ratio_diff_mother.Sumw2()
    k_star_ratio_diff_mother.Rebin(10)
    k_star_ratio_diff_mother.GetXaxis().SetTitle("#it{k*} (GeV/#it{c})")
    k_star_ratio_diff_mother.GetYaxis().SetTitle("C(#it{k*})")

    thnsparse_same.GetAxis(3).SetRange(thnsparse_same.GetAxis(3).FindBin(0.5), thnsparse_same.GetAxis(3).FindBin(0.5))
    k_star_ratio_same_mother = thnsparse_same.Projection(4).Clone(f"k_star_ratio_same_mother_{baryon_type}")
    k_star_ratio_same_mother.Sumw2()
    k_star_ratio_same_mother.Rebin(10)
    k_star_ratio_same_mother.GetXaxis().SetTitle("#it{k*} (GeV/#it{c})")
    k_star_ratio_same_mother.GetYaxis().SetTitle("C(#it{k*})")

    mix_events = thnsparse_mix.Projection(4)
    mix_events.SetName(f"mix_events_{baryon_type}")
    mix_events.Sumw2()
    mix_events.Rebin(10)
    mix_events.Scale(1.0/n_mix)
    k_star_ratio_inclusive.Divide(mix_events)
    k_star_ratio_diff_mother.Divide(mix_events)
    k_star_ratio_same_mother.Divide(mix_events)
    k_star_ratio_inclusive.Write()
    k_star_ratio_diff_mother.Write()
    k_star_ratio_same_mother.Write()

out_file.cd()
h_sigma_pt.Write("sigma_pt_distribution")
h_proton_pt.Write("proton_pt_distribution")
out_file.Close()