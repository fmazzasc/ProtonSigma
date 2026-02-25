import ROOT
ROOT.gROOT.SetBatch(True)

def group_bins_fixed(ny, group_size):
    groups = []
    start = 1
    while start <= ny:
        end = min(ny, start + group_size - 1)
        groups.append(list(range(start, end + 1)))
        start = end + 1
    return groups

def x_center_for_group(yaxis, ybins, mode="mean"):
    centers = [yaxis.GetBinCenter(b) for b in ybins]
    if mode == "mid":
        return 0.5 * (centers[0] + centers[-1])
    return sum(centers) / len(centers)

def project_y_sum_xbins(h2, name, xbins):
    hsum = None
    for xb in xbins:
        tmp = h2.ProjectionY(f"{name}_xb{xb}", xb, xb)
        tmp.SetDirectory(0)
        if hsum is None:
            hsum = tmp.Clone(name)
            hsum.SetDirectory(0)
        else:
            hsum.Add(tmp)
    return hsum

file_in = ROOT.TFile.Open("kstar_sigma_proton_analysis_MC_recal.root")
h2 = file_in.Get("h2KStarResolutionVsKstar")
h2.SetDirectory(0)

file_data = ROOT.TFile.Open("histo_subtracted.root")
me_histo = file_data.Get("mixed_event")
me_histo.SetDirectory(0)


output_file = ROOT.TFile.Open("smeared_CF.root", "RECREATE")
print('-------------------------------')
print("Get resolution vs k* from MC, Gausssian fit of the slices")
## make groouped slices
groups = group_bins_fixed(h2.GetNbinsX(), group_size=5)
mean_arr = []
mean_err_arr = []
sigma_arr = []
sigma_err_arr = []
k_center_arr = []
outdir = output_file.mkdir("slices")
outdir.cd()
for ig, ybins in enumerate(groups, start=1):
    k_center_arr.append(x_center_for_group(h2.GetXaxis(), ybins, mode="mean"))
    proj = project_y_sum_xbins(h2, f"proj_{ig}", ybins)
    if proj.GetEntries() < 10:
        mean_arr.append(0)
        mean_err_arr.append(0)
        sigma_arr.append(0)
        sigma_err_arr.append(0)
        continue

    fit_result = proj.Fit("gaus", "QRMS", "", -0.1, 0.1)
    proj.Write()
    mean = fit_result.Parameter(1)
    mean_err = fit_result.ParError(1)
    sigma = fit_result.Parameter(2)
    sigma_err = fit_result.ParError(2)

    mean_arr.append(mean)
    mean_err_arr.append(mean_err)
    sigma_arr.append(sigma)
    sigma_err_arr.append(sigma_err)


## create output histograms for mean and sigma vs k*
hMean = ROOT.TH1D("hMean", "Mean vs k*; k* (GeV/c); Mean (GeV/c)", len(groups), 0, max(k_center_arr)*1.1)
hSigma = ROOT.TH1D("hSigma", "Sigma vs k*; k* (GeV/c); Sigma (GeV/c)", len(groups), 0, max(k_center_arr)*1.1)
for i, (k, m, me, s, se) in enumerate(zip(k_center_arr, mean_arr, mean_err_arr, sigma_arr, sigma_err_arr), start=1):
    hMean.SetBinContent(i, m)
    hMean.SetBinError(i, me)
    hSigma.SetBinContent(i, s)
    hSigma.SetBinError(i, se)
## convert gev (res) in mev (cf)
hMean.Scale(1000)
hSigma.Scale(1000)
## fit mean and sigma as a function of k* with a polynomial
fit_func = ROOT.TF1("pol3", "pol3", 0, max(k_center_arr)*1.1)
hMean.Fit(fit_func, "QRMS", "", 0, 1)
pars_mean = [hMean.GetFunction("pol3").GetParameter(i) for i in range(4)]
hSigma.Fit(fit_func, "QRMS", "", 0, 1)
pars_sigma = [hSigma.GetFunction("pol3").GetParameter(i) for i in range(4)]
## get parameters of the fits
print('-------------------------------')
print("Fit mean and sigma vs k* with a polynomial, get fit parameters")
print("Mean fit parameters:", pars_mean)
print("Sigma fit parameters:", pars_sigma)

## start smearing the CF by convolving with a Gaussian whose mean and sigma depend on k* according to the fitted functions
output_file.cd()
hMean.Write("hMean_fit")
hSigma.Write("hSigma_fit")

theo_file = ROOT.TFile.Open("theory_CFs_all.root")
theo_dirs = []
for key in theo_file.GetListOfKeys():
    if key.IsFolder():
        theo_dirs.append(key.GetName())

for dir in theo_dirs:
    print("--------------------------------")
    print(f"Smearing CF from: {dir}")
    corr = theo_file.Get(f"{dir}/CF_Combined")
    n_gr_points = corr.GetN()
    smeared_graph = ROOT.TGraph(n_gr_points)
    for i in range(n_gr_points):
        k = corr.GetX()[i]
        gaus = ROOT.TF1("gaus", "gaus", -2, 2)
        gaus_mean = pars_mean[0] + pars_mean[1] * k / 1000 + pars_mean[2] * (k / 1000)**2 + pars_mean[3] * (k / 1000)**3 + k 
        gaus_sigma = pars_sigma[0] + pars_sigma[1] * k / 1000 + pars_sigma[2] * (k / 1000)**2 + pars_sigma[3] * (k / 1000)**3
        gaus.SetParameters(1, gaus_mean, gaus_sigma)
        # Perform the convolution integral for the CF at this k
        integral = 0
        norm = 0
        for j in range(n_gr_points):
            k_prime = corr.GetX()[j]
            cf_value = corr.GetY()[j]
            me_value = me_histo.GetBinContent(me_histo.FindBin(k_prime / 1000))
            weight = gaus.Eval(k_prime) * me_value
            integral += cf_value * weight
            norm += weight
        smeared_value = integral / norm if norm > 0 else 1.0
        smeared_graph.SetPoint(i, k, smeared_value)
    outdir = output_file.mkdir(dir + "_smeared")
    outdir.cd()
    corr.Write("original_CF")
    smeared_graph.Write("smeared_CF")