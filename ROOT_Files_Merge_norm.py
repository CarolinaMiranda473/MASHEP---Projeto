import ROOT
import os
import glob

##This codes joins all the Monte Carlo Files that correspond to the background of ZZ  onto the same file (ZZBackground.root),
#and joind sub-leading background onto the Root file MinorBackgrounds.root.
##It also joins all the signal files onto the same file (Signal.root), regardless of Higgs production mode, since all will be analysed 
# simultaneously.
##During it, we make sure that the MC histograms are normalzied to the number of expected events for each process.

# List of histogram names - All files contain the same histograms
histogram_names = [
    "hist_mLL1",
    "hist_mLL2",
    "hist_fourlepsys_pt",
    "hist_fourlepsys_y",
    "mass_four_lep",
    "mass_ext_four_lep",
    "hist_n_jets",
    "hist_fourleptpt",
    "hist_fourlepteta",
    "hist_fourleptE",
    "hist_fourleptphi",
    "hist_fourleptID"
]

lumi = 10064 ##pb-1

# Normalization factors for each histogram
#Determined by Norm_factor = N_expected/N_simulated, where 
# N_expected = cross section * luminosity; N_simulated = sum of weights [sum of weights and cross section obtained from Files_HZZ.txt]
#Carfeul: Z_Z->Scale(1.3); // loop-induced gluon–gluon gg->ZZ is not included in the current MCs, llll, ZqqZll, llvv
normalization_factors = {
    "Zee.root": 1950.5295 * lumi / 150277594200,
    "Zmumu.root": 1950.6321 * lumi / 147334691090,
    "Ztautau.root": 1950.6321 * lumi / 56171652547.3,
    "Wplusenu.root": 11500.4632 * lumi / 473389396815,
    "Wplusmunu.root": 11500.4632 * lumi / 446507925520,
    "Wplustaunu.root": 11500.4632 * lumi / 670928468875,
    "Wminusenu.root": 8579.63498 * lumi / 247538642447,
    "Wminusmunu.root": 8579.63498 * lumi / 264338188182,
    "Wminustaunu.root": 8579.63498 * lumi / 165195850954,
    "ZqqZll.root": 2.20355112 * lumi * 1.3/ 3439266.11559,
    "WqqZll.root": 3.4328 * lumi / 241438.72705,
    "WpqqWmlv.root": 24.708 * lumi / 998250.783475,
    "WplvWmqq.root": 24.724 * lumi / 1069526.41899,
    "WlvZqq.root": 11.42 * lumi / 1111991.15979,
    "llll.root": 1.2578 * lumi * 1.3 / 7538705.8077,
    "lllv.root": 4.6049 * lumi / 5441475.00407,
    "llvv.root": 12.466 * lumi * 1.3 / 5039259.9696,
    "lvvv.root": 3.2286 * lumi / 1727991.07441,
    "ttbar_lep.root": 452.693559 * lumi / 49386600,
    "single_top_tchan.root": 44.152 * lumi / 4986200,
    "single_antitop_tchan.root": 26.276 * lumi / 4989800,
    "single_top_schan.root": 2.06121 * lumi / 997800,
    "single_antitop_schan.root": 1.288662 * lumi / 995400,
    "single_top_wtchan.root": 35.845486 * lumi / 4865800,
    "single_antitop_wtchan.root": 35.824406 * lumi / 4945600,
    "ggH125_ZZ4lep.root": 0.0060239 * lumi / 27881776.6536,
    "ZH125_ZZ4lep.root": 0.0000021424784 * lumi / 150000,
    "WH125_ZZ4lep.root": 0.0003769 * lumi / 149400,
    "VBFH125_ZZ4lep.root": 0.0004633012 * lumi / 3680490.83243,
}
# Directory containing ROOT files
input_dir = "Output_HZZAnalysis"
root_files = glob.glob(os.path.join(input_dir, "*.root"))

# Exclude specific files (e.g., data)
excluded_files = {"data.root"}

# Identify signal files
signal_files = [f for f in root_files if "H125" in f and os.path.basename(f) not in excluded_files]

# Remaining files are considered backgrounds
background_files = [
    f for f in root_files
    if f not in signal_files  and os.path.basename(f) not in excluded_files
]

def Merge_Histograms_Func(file_list, output_filename):
    summed_hists = {name: None for name in histogram_names}

    for file in file_list:
        file_basename = os.path.basename(file)
        factor = normalization_factors.get(file_basename, 1.0)  # default to 1.0 if not specified

        f = ROOT.TFile.Open(file)
        if not f or f.IsZombie():
            print(f"Could not open {file}, skipping.")
            continue

        for name in histogram_names:
            hist = f.Get(name)
            if not hist:
                print(f"Histogram {name} not found in {file}, skipping.")
                continue

            hist.SetDirectory(0)
            hist.Scale(factor)

            if summed_hists[name] is None:
                summed_hists[name] = hist.Clone()
                summed_hists[name].SetDirectory(0)
            else:
                summed_hists[name].Add(hist)

        f.Close()

    output_path = os.path.join(input_dir, output_filename)
    out_file = ROOT.TFile(output_path, "RECREATE")
    for hist in summed_hists.values():
        if hist:
            hist.Write()
    out_file.Close()
    print(f"Written: {output_path}")

# Merge and write signal and background ROOT files
Merge_Histograms_Func(signal_files, "Signal.root")
Merge_Histograms_Func(background_files, "Background.root")
