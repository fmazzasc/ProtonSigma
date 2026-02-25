import ROOT
from array import array

ROOT.gROOT.SetBatch(True)

# -----------------------
# Helpers
# -----------------------
def integral_in_range(h, x_min, x_max):
    b1 = h.FindBin(x_min)
    b2 = h.FindBin(x_max)
    return h.Integral(b1, b2)

def integral_error_in_range(h, x_min, x_max):
    b1 = h.FindBin(x_min)
    b2 = h.FindBin(x_max)
    error2 = sum(h.GetBinError(i)**2 for i in range(b1, b2 + 1))
    return error2 ** 0.5

def project_x_sum_ybins(h2, name, ybins):
    """
    Returns a TH1 projection on X, summing over the provided list of Y-bin indices.
    """
    hsum = None
    for yb in ybins:
        tmp = h2.ProjectionX(f"{name}_yb{yb}", yb, yb)
        tmp.SetDirectory(0)
        if hsum is None:
            hsum = tmp.Clone(name)
            hsum.SetDirectory(0)
        else:
            hsum.Add(tmp)
    return hsum

def group_bins_fixed(ny, group_size):
    """
    Groups Y bins in contiguous chunks of size group_size.
    group_size=1 reproduces the original behavior.
    """
    groups = []
    start = 1
    while start <= ny:
        end = min(ny, start + group_size - 1)
        groups.append(list(range(start, end + 1)))
        start = end + 1
    return groups

def group_bins_custom(list_of_groups):
    """
    User-provided groups like [[1,2,3],[4],[5,6]].
    """
    return list_of_groups

def y_center_for_group(yaxis, ybins, mode="mean"):
    centers = [yaxis.GetBinCenter(b) for b in ybins]
    if mode == "mid":
        return 0.5 * (centers[0] + centers[-1])
    return sum(centers) / len(centers)

def subtract_background(
    proj_data, proj_me, proj_mc_bkg,
    norm_min=1.24, norm_max=1.30,
    comp_dir=None, data_dir=None, me_dir=None, bkg_dir=None,
    tag=""
):
    """
    Normalize MC bkg to SE in [norm_min, norm_max], save bkg projections,
    subtract from SE, rescale bkg to ME level, save again, subtract from ME.
    """
    data_int = integral_in_range(proj_data, norm_min, norm_max)
    me_int   = integral_in_range(proj_me,   norm_min, norm_max)
    mc_int   = integral_in_range(proj_mc_bkg, norm_min, norm_max)

    bin_width_data = proj_data.GetXaxis().GetBinWidth(1)
    bin_width_me   = proj_me.GetXaxis().GetBinWidth(1)
    bin_width_mc   = proj_mc_bkg.GetXaxis().GetBinWidth(1)

    if data_int <= 0 or me_int <= 0 or mc_int <= 0:
        return None, None

    # Normalize MC background to Same Event
    proj_mc_bkg.Scale(data_int / mc_int * (bin_width_data / bin_width_mc))

    # Save background as used for SE subtraction
    if bkg_dir is not None:
        bkg_dir.cd()
        proj_mc_bkg.Write(f"bkg_norm_SE_{tag}")

    # Comparison canvas (SE vs bkg)
    if comp_dir is not None:
        cv = ROOT.TCanvas(f"cv_{tag}", f"Comparison {tag}", 800, 600)
        proj_data.SetLineColor(ROOT.kBlue)
        proj_data.Draw("E")
        proj_mc_bkg.SetLineColor(ROOT.kRed)
        proj_mc_bkg.Draw("E SAME")
        comp_dir.cd()
        cv.Write()

    # Subtract from Same Event
    proj_data.Add(proj_mc_bkg, -1)
    if data_dir is not None:
        data_dir.cd()
        proj_data.Write(f"proj_data_sub_{tag}")

    mc_int_after = integral_in_range(proj_mc_bkg, norm_min, norm_max)
    proj_mc_bkg.Scale(me_int / mc_int_after * (bin_width_me / bin_width_mc))

    # Subtract from Mixed Event
    proj_me.Add(proj_mc_bkg, -1)
    if me_dir is not None:
        me_dir.cd()
        proj_me.Write(f"proj_me_sub_{tag}")
    
    # Save background as used for ME subtraction
    if bkg_dir is not None:
        bkg_dir.cd()
        ## normalize to unity background
        proj_mc_bkg.Scale(1.0 / proj_mc_bkg.Integral())
        proj_mc_bkg.Write(f"bkg_norm_ME_{tag}")

    return proj_data, proj_me

def make_variable_edges_from_centers(centers):
    """
    Build variable bin edges from bin centers (sorted).
    Edges are midpoints between adjacent centers, with outer edges extrapolated.
    """
    c = sorted(centers)
    if len(c) == 1:
        return [c[0] - 0.5, c[0] + 0.5]

    mids = [(c[i] + c[i+1]) * 0.5 for i in range(len(c) - 1)]
    left_width  = mids[0] - c[0]
    right_width = c[-1] - mids[-1]
    edges = [c[0] - left_width] + mids + [c[-1] + right_width]
    return edges

# -----------------------
# Inputs
# -----------------------
file_data = ROOT.TFile("kstar_sigma_proton_analysis.root")
histo_data = file_data.Get("hKStarVsMassSigma")
histo_data.SetDirectory(0)

file_data_me = ROOT.TFile("kstar_sigma_proton_analysis_ME.root")
histo_data_me = file_data_me.Get("hKStarVsMassSigma")
histo_data_me.SetDirectory(0)

file_mc = ROOT.TFile("kstar_sigma_proton_analysis_MC.root")
histo_mc_bkg = file_mc.Get("h2KstarMassSigmaBkg")
histo_mc_bkg.SetDirectory(0)

# -----------------------
# Output structure
# -----------------------
outfile = ROOT.TFile("histo_subtracted.root", "RECREATE")
comp_dir = outfile.mkdir("comparisons")
data_dir = outfile.mkdir("projections_data")
me_dir   = outfile.mkdir("projections_me")
bkg_dir  = outfile.mkdir("projections_background")

# -----------------------
# Choose grouping of Y bins
# -----------------------
ny = histo_data.GetNbinsY()

# Option A: original behavior
y_groups = group_bins_fixed(ny, group_size=3)
k_star_centers = []
same_event_counts = []
same_event_errors = []
mixed_event_counts = []
mixed_event_errors = []

print("histo MC entries:", histo_mc_bkg.GetEntries())
bin_y_08 = histo_mc_bkg.GetYaxis().FindBin(0.9)
proj_mc_08 = project_x_sum_ybins(histo_mc_bkg, "proj_mc_08", list(range(1, bin_y_08 + 1)))
bin_y_13 = histo_mc_bkg.GetYaxis().FindBin(1.3)
proj_mc_13 = project_x_sum_ybins(histo_mc_bkg, "proj_mc_13", list(range(bin_y_08 + 1, bin_y_13 + 1)))
bin_y_max = histo_mc_bkg.GetNbinsY()
proj_mc_13p = project_x_sum_ybins(histo_mc_bkg, "proj_mc_13p", list(range(bin_y_13 + 1, bin_y_max + 1)))

proj_mc_bkg = project_x_sum_ybins(histo_mc_bkg, "proj_mc_bkg_all", list(range(1, ny + 1)))
for ig, ybins in enumerate(y_groups, start=1):
    tag = f"grp{ig}_y{'_'.join(map(str, ybins))}"
    kstar_bin_center = y_center_for_group(histo_data.GetYaxis(), ybins)
    print(f"Group {ig}: K* bin center = {kstar_bin_center:.4f}, MC bkg integral in norm range = {integral_in_range(proj_mc_bkg, 1.27, 1.3):.0f}")

    proj_data = project_x_sum_ybins(histo_data,   f"proj_data_{tag}", ybins)
    proj_me = project_x_sum_ybins(histo_data_me,f"proj_me_{tag}",  ybins)
    if proj_mc_bkg is None or proj_data is None or proj_me is None:
        continue
    
    # print(f"Group {ig}: K* bin center = {kstar_bin_center:.4f}, MC bkg integral in norm range = {integral_in_range(proj_mc_bkg, 1.27, 1.3):.0f}")

    proj_data_sub, proj_me_sub = subtract_background(
        proj_data, proj_me, proj_mc_bkg,
        norm_min=1.27, norm_max=1.3,
        comp_dir=comp_dir,
        data_dir=data_dir,
        me_dir=me_dir,
        bkg_dir=bkg_dir,
        tag=tag
    )
    if proj_data_sub is None:
        continue

    kstar_c = y_center_for_group(histo_data.GetYaxis(), ybins, mode="mean")
    k_star_centers.append(kstar_c)
    se = integral_in_range(proj_data_sub, 1.18, 1.22)
    me = integral_in_range(proj_me_sub,   1.18, 1.22)
    se_error = integral_error_in_range(proj_data_sub, 18, 1.22)
    me_error = integral_error_in_range(proj_me_sub,   1.18, 1.22)
    print(f"Group {ig}: K* center = {kstar_c:.4f}, SE = {se:.0f} ± {se_error:.0f}, ME = {me:.0f} ± {me_error:.0f}")
    same_event_counts.append(se)
    mixed_event_counts.append(me)
    same_event_errors.append(se_error)
    mixed_event_errors.append(me_error)

# -----------------------
# Ratio histogram
# -----------------------
if len(k_star_centers) > 0:
    # Build variable bin edges from centers (robust for grouped bins)
    edges = make_variable_edges_from_centers(k_star_centers)
    arr_edges = array('d', edges)

    ratio_hist = ROOT.TH1D("ratio", "Same Event / Mixed Event Ratio;K*;SE/ME", len(arr_edges) - 1, arr_edges)
    se_hist = ROOT.TH1D("same_event", "Same Event Counts;K*;Counts", len(arr_edges) - 1, arr_edges)
    me_hist = ROOT.TH1D("mixed_event", "Mixed Event Counts;K*;Counts", len(arr_edges) - 1, arr_edges)

    # Fill bins in order of increasing center
    order = sorted(range(len(k_star_centers)), key=lambda i: k_star_centers[i])
    for ibin, idx in enumerate(order, start=1):
        c  = k_star_centers[idx]
        se = same_event_counts[idx]
        me = mixed_event_counts[idx]
        se_error = same_event_errors[idx]
        me_error = mixed_event_errors[idx]
        se_hist.SetBinContent(ibin, se)
        se_hist.SetBinError(ibin, se_error)
        me_hist.SetBinContent(ibin, me)
        me_hist.SetBinError(ibin, me_error)

        if se > 0 and me > 0:
            ratio_hist.SetBinContent(ibin, se / me)
            error = (se_error / se)**2 + (me_error / me)**2
            ratio_hist.SetBinError(ibin, (se / me) * error**0.5)
        else:
            ratio_hist.SetBinContent(ibin, 0)
            ratio_hist.SetBinError(ibin, 0)
        print(f"K*: {c:.4f}, Same Event: {se:.0f}, Mixed Event: {me:.0f}, Ratio: {ratio_hist.GetBinContent(ibin):.4f} ± {ratio_hist.GetBinError(ibin):.4f}")

    outfile.cd()
    me_hist.Write()
    se_hist.Write()
    ratio_hist.Write()

outfile.Close()
