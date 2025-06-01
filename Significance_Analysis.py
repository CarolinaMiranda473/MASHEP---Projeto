import ROOT
import math
import numpy as np
from ROOT.Math import poisson_pdf
import scipy


##First, we shall define some basic functions, that will be used throughout our analysis.

## ------------------------------FUNCTION DEFINITIONS-------------------------------------------------
#Function to calculate Significance - For counting Statistics - More accurate
def Significance_Func_1(sig, bkg):
    if bkg == 0:
        return 0
    else:
        Sig_result = math.sqrt(2 * (sig + bkg) * math.log(1 + sig/bkg) - 2 * sig)
    
    return Sig_result

#Function to calculate Significance - For counting Statistics - Less accurate
def Significance_Func_2(sig, bkg):
    
    Sig_result = sig / math.sqrt(sig + bkg)
    
    return Sig_result

#Function to calculate selection Efficiency 
def Efficiency_Func(Ntot, Nsel):

    if Ntot <= 0:
        return 0
    
    else:
        return Nsel / Ntot 
    
#Function to calculate Log Likelihood
def Log_Likelihood_Func(mu, signal_exp, bkg_exp, data):
  
    likelihood_result = np.zeros(signal_exp.GetNbinsX())
    nbins = signal_exp.GetNbinsX()

    for i in range(nbins):
        expected = mu*signal_exp.GetBinContent(i) + bkg_exp.GetBinContent(i)
        datavalue = data.GetBinContent(i)
        if expected != 0:
            likelihood_result[i] = math.log(poisson_pdf(int(datavalue),expected)) 
        else:
            likelihood_result[i] = 0

    total_ll = np.sum(likelihood_result)

    return total_ll

def compute_sigma_intervals(mu_values, log_likelihood):
    import numpy as np

    # Convert to numpy arrays
    mu_values = np.array(mu_values)
    log_likelihood = np.array(log_likelihood)

    # Get the minimum value and subtract to find ΔlogL
    minNLogValue = np.max(log_likelihood)  # Note: log-likelihood is usually negative
    Delta_neglogL = minNLogValue - log_likelihood

    # Confidence level thresholds (ΔlogL)
    CL_thresholds = {
        "1σ": 0.5,
        "2σ": 2.0,
        "3σ": 4.5
    }

    intervals = {}

    for label, threshold in CL_thresholds.items():
        lower_mu = upper_mu = None

        for i in range(1, len(mu_values)):
            # Look for crossing from below to above (left to right)
            if Delta_neglogL[i-1] < threshold <= Delta_neglogL[i]:
                x0, x1 = mu_values[i-1], mu_values[i]
                y0, y1 = Delta_neglogL[i-1], Delta_neglogL[i]
                lower_mu = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)

            # Look for crossing from above to below (right to left)
            if Delta_neglogL[i-1] > threshold >= Delta_neglogL[i]:
                x0, x1 = mu_values[i-1], mu_values[i]
                y0, y1 = Delta_neglogL[i-1], Delta_neglogL[i]
                upper_mu = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)

        if lower_mu is not None and upper_mu is not None:
            if upper_mu < lower_mu:
                lower_mu, upper_mu = upper_mu, lower_mu

            intervals[label] = {
                "lower": lower_mu,
                "upper": upper_mu,
                "half_width": (upper_mu - lower_mu) / 2
            }
        else:
            intervals[label] = {
                "lower": None,
                "upper": None,
                "half_width": None
            }

    return intervals

## --------------------------------IMPORTING FILES---------------------------------------------------
##Reading all the files
#We have joined all background types into only one Background.root file, with all the needed histograms.
#We have also joined all Higgs production types' Signals into only one Signal.root file, with all the needed histograms.
#The real data is in the data.root file, taken directly from the ATLAS open data framework.~

# Open the ROOT files
Signal_file = ROOT.TFile("./Output_HZZAnalysis/Signal.root", "READ")
Background_file = ROOT.TFile("./Output_HZZAnalysis/Background.root", "READ")
Data_file = ROOT.TFile("./Output_HZZAnalysis/data.root", "READ")

# List of histograms contained in each file
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

# List of variable names, to be used in the X-axis of our histograms
x_axis_titles = [
    "m_{ll}^{(1)} [GeV]",         # hist_mLL1
    "m_{ll}^{(2)} [GeV]",         # hist_mLL2
    "p_{T}^{4l system} [GeV]",    # hist_fourlepsys_pt
    "Rapidity of 4l system",     # hist_fourlepsys_y
    "m_{4l} [GeV]",              # mass_four_lep
    "Extended m_{4l} [GeV]",     # mass_ext_four_lep
    "Jet multiplicity",          # hist_n_jets
    "p_{T}^{4l} [GeV]",          # hist_fourleptpt
    "η^{4l}",                    # hist_fourlepteta
    "E^{4l} [GeV]",              # hist_fourleptE
    "ϕ^{4l}",                    # hist_fourleptphi
    "4l PDG ID sum"              # hist_fourleptID
]

# List of histogram titles
plot_titles = [
    "First Lepton Pair Invariant Mass",         # hist_mLL1
    "Second Lepton Pair Invariant Mass",        # hist_mLL2
    "Transverse Momentum of 4-lepton system",   # hist_fourlepsys_pt
    "Rapidity of 4-lepton system ",             # hist_fourlepsys_y
    "4-lepton Invariant Mass",                  # mass_four_lep
    "Extended 4-lepton Invariant Mass",         # mass_ext_four_lep
    "Jet multiplicity",                         # hist_n_jets
    "Transverse Momentum of 4-lepton",          # hist_fourleptpt
    "Rapidity of 4-lepton",                     # hist_fourlepteta
    "Energy of 4-lepton",                       # hist_fourleptE
    "Phi of 4-lepton",                          # hist_fourleptphi
    "4l PDG ID sum"                             # hist_fourleptID
]

hists_len = len(histogram_names)

## -----------------------------------DRAWING HISTOGRAMS--------------------------------------------------
c = ROOT.TCanvas()

# Initialize lists to hold histograms and legends
sig_hist = [None] * hists_len
bkg_hist = [None] * hists_len
data_hist = [None] * hists_len
legends = []

# Function to draw Histograms -> Maybe change colors to be different for each observable, and not always the same
# Used TStack to sum Signal + Background, to be plotted against data -> Same as used in the HZZAnalysis folder!
def Drawing(i, histo_name, x_axis_title="", plot_title=""):
    # Retrieve histograms
    sig_hist[i] = Signal_file.Get(histo_name)
    bkg_hist[i] = Background_file.Get(histo_name)
    data_hist[i] = Data_file.Get(histo_name)

    mc_stack = ROOT.THStack("mc_stack", f"{plot_title}; {x_axis_title}; Events")

    # Check for missing histograms
    if not sig_hist[i] or not bkg_hist[i] or not data_hist[i]:
        print(f"Error: Histogram '{histo_name}' not found in one of the files.")
        return

    # Histogram aesthetics - colors and styles
    sig_hist[i].SetLineColor(ROOT.kRed-4)
    sig_hist[i].SetFillColor(ROOT.kRed-4)
    sig_hist[i].SetFillStyle(3001)

    bkg_hist[i].SetLineColor(ROOT.kAzure+8)
    bkg_hist[i].SetFillColor(ROOT.kAzure+8)
    bkg_hist[i].SetFillStyle(3002)

    data_hist[i].SetLineColor(ROOT.kBlack)
    data_hist[i].SetMarkerStyle(20)
    data_hist[i].SetMarkerColor(ROOT.kBlack)
    data_hist[i].SetLineWidth(2)

    sig_hist[i].SetStats(0)
    bkg_hist[i].SetStats(0)
    data_hist[i].SetStats(0)

    # Get max value of data histogram and scale it
    ymax = data_hist[i].GetMaximum() * 1.3

    # Draw as histograms
    #bkg_hist[i].Draw("hist")
    #sig_hist[i].Draw("histsame")
    #data_hist[i].Draw("esame")


    #Draw using a TStack -> Sums Signal and Background contributions
    mc_stack.SetMaximum(ymax)
    mc_stack.Add(bkg_hist[i])
    mc_stack.Add(sig_hist[i])
    mc_stack.Draw("hist")
    data_hist[i].Draw("esame")
    

    # Legend
    leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    leg.SetHeader(histo_name, "C")
    leg.AddEntry(sig_hist[i], "Signal", "f")
    leg.AddEntry(bkg_hist[i], "Backgrounds", "f")
    leg.AddEntry(data_hist[i], "Data", "lep")
    leg.Draw()

    output_path = f"Results_Histos/{histo_name}_Ours.pdf"
    
    c.Print(output_path)

    c.Clear()


#Draw all the acessible variables
for i in range(hists_len):
    Drawing(
        i,
        histogram_names[i],
        x_axis_titles[i],
        plot_titles[i]
    )

## -----------------------------------------CHOOSING SELECTION VARIABLES ------------------------------------------
#We choose the variable that best optimizes signal significance, and apply cuts on it, to perform a discovery fit.

best_index = list(range(hists_len))

#Create list to store Significance histograms 

Significance_Right_Cut = []
Significance_Left_Cut = []

c.Divide(6, 2)


##Determining best Signal Cut for each variable available

for i, var_index in enumerate(best_index):
    # Clone the histograms
    sig_hist_clone = sig_hist[var_index].Clone()
    bkg_hist_clone = bkg_hist[var_index].Clone()

    # Right cut significance
    significance_right = sig_hist_clone.Clone(f"{histogram_names[var_index]}_right")
    significance_right.Reset()
    n_bins = sig_hist_clone.GetNbinsX()

    for j in range(1, n_bins + 1):
        S_right = sig_hist_clone.Integral(1, j)
        B_right = bkg_hist_clone.Integral(1, j)
        significance_right.SetBinContent(j, Significance_Func_1(S_right, B_right))

    Significance_Right_Cut.append(significance_right)

    # Left cut significance
    significance_left = sig_hist_clone.Clone(f"{histogram_names[var_index]}_left")
    significance_left.Reset()
    
    for j in range(1, n_bins + 1):
        S_left = sig_hist_clone.Integral(j, n_bins)
        B_left = bkg_hist_clone.Integral(j, n_bins)
        significance_left.SetBinContent(j, Significance_Func_1(S_left, B_left))

    Significance_Left_Cut.append(significance_left)

    # Print initial total significance
    initial_significance_right = Significance_Right_Cut[i].Integral()
    initial_significance_left = Significance_Left_Cut[i].Integral()
    print(f"Initial Total Significance for {histogram_names[var_index]}:")
    print(f"  Right Cut: {initial_significance_right:.3f}")
    print(f"  Left Cut: {initial_significance_left:.3f}")

colors = [ ##Actually need to choose colors -> Different color for each variable: Could also be done for the first Histograms!!THINK ABOUT IT
    ROOT.kRed + 1, 
    ROOT.kRed + 2, 
    ROOT.kRed + 3, 
    ROOT.kRed + 4, 
    ROOT.kRed + 5, 
    ROOT.kRed + 6, 
    ROOT.kRed + 7,
    ROOT.kRed + 8, 
    ROOT.kRed + 9,
    ROOT.kRed + 10,
    ROOT.kPink +1, 
    ROOT.kPink + 2
    ]

## Plotting Significance value vs Cut on variable, for all available values
# Draw Right Cut Significance
for i, sig in enumerate(Significance_Right_Cut):
    c.cd(i + 1)
    sig.SetLineColor(colors[i])
    sig.SetFillColor(ROOT.kGray)
    sig.SetFillStyle(3001)
    sig.SetStats(0)
    sig.GetXaxis().SetTitle(f"{x_axis_titles[best_index[i]]} [GeV]")
    sig.GetYaxis().SetTitle(f"Significance 0 -> i")
    sig.Draw("hist")

c.Print("SignificanceRight_0toi.pdf")
c.Clear()


# Draw Left Cut Significance
c.Divide(6, 2)
for i, sig in enumerate(Significance_Left_Cut):
    c.cd(i + 1)
    sig.SetLineColor(colors[i])
    sig.SetFillColor(ROOT.kGray)
    sig.SetFillStyle(3001)
    sig.SetStats(0)
    sig.GetXaxis().SetTitle(f"{x_axis_titles[best_index[i]]} [GeV]")
    sig.GetYaxis().SetTitle(f"Significance i -> N")
    sig.Draw("hist")

c.Print("SignificanceLeft_itoN.pdf")
c.Clear()

# Print results
left_cuts = [0.0] * len(best_index)
right_cuts = [0.0] * len(best_index)

#Now, we look at the maximum significances achieved, with Left to Right and Right to Left cuts, for each variable.
# If we only want to do univariate optimization (that is, only use one variable, and apply the best cut on it) we must only look
# at the maximum value of significance, to identify: for which variable does it occur; for what cut values.
print(f"--------------------------------------BEST CUTS---------------------------------------")
for i, var_index in enumerate(best_index):
    left_cuts[i] = Significance_Left_Cut[i].GetXaxis().GetBinCenter(Significance_Left_Cut[i].GetMaximumBin())
    right_cuts[i] = Significance_Right_Cut[i].GetXaxis().GetBinCenter(Significance_Right_Cut[i].GetMaximumBin())
    print(f"{histogram_names[var_index]}: {left_cuts[i]:.3f} -> {right_cuts[i]:.3f}")

    # Print final total significance after cuts
    sighist_clone = sig_hist[var_index]
    bkghist_clone = bkg_hist[var_index]
    left_bin = sighist_clone.GetXaxis().FindBin(left_cuts[i])
    right_bin = sighist_clone.GetXaxis().FindBin(right_cuts[i])

    S_left = sighist_clone.Integral(left_bin, sighist_clone.GetNbinsX())
    B_left = bkghist_clone.Integral(left_bin, bkghist_clone.GetNbinsX())

    S_right = sighist_clone.Integral(1, right_bin)
    B_right = bkghist_clone.Integral(1, right_bin)

    final_significance_left1 = Significance_Func_1(S_left, B_left)
    final_significance_right1 = Significance_Func_1(S_right, B_right)

    final_significance_left2 = Significance_Func_2(S_left, B_left)
    final_significance_right2 = Significance_Func_2(S_right, B_right)
    print(f"Final Total Significance for {histogram_names[best_index[i]]}:")
    print(f"  Left Cut: Function 1 (More Accurate) -> {final_significance_left1:.3f}")
    print(f"  Left Cut: Function 2 -> {final_significance_left2:.3f}")
    print(f"  Right Cut: Function 1 (More Accurate) -> {final_significance_right1:.3f}")
    print(f"  Right Cut: Function 2 -> {final_significance_right2:.3f}")

## -----------------------------------------PERFORMING SIGNAL SELECTION CUTS---------------------------------------
#We perform signal selection, based on maximum significance values taken from Monte Carlo simulations, to optimize our max 
#likelihood fits (We can also use ML to optimize this fit -> We will do so, in a different analysis strategy, in NOME.py .)

#By looking at the values obtained, we find that:
# Optimal variable to apply cut: mass_four_lep
# Optimal cut interval: 111.167 -> 128.167

#Name and index of variable we will perform selection cuts on.
opt_var = "mass_four_lep" 
opt_var_index = 4


# Creating histograms, after signal selection. We will iterate over the histograms, and make every bin outside our selection 0.
sig_hist_SEL = sig_hist[opt_var_index].Clone(f"{opt_var}_sigSEL")
bkg_hist_SEL = bkg_hist[opt_var_index].Clone(f"{opt_var}_bkgSEL")
data_hist_SEL = data_hist[opt_var_index].Clone(f"{opt_var}_dataSEL")

def Cut_Window(hist, left_cut, right_cut):

    for bin in range(1, hist.GetNbinsX() + 1):
        bin_center = hist.GetXaxis().GetBinCenter(bin)
        
        if bin_center < left_cut or bin_center > right_cut:
            hist.SetBinContent(bin, 0)

for hist in [sig_hist_SEL, bkg_hist_SEL, data_hist_SEL]:
    Cut_Window(hist, left_cuts[opt_var_index], right_cuts[opt_var_index])

## Plotting Histograms for 4-lepton invariant mass BEFORE vs AFTER signal selection cuts

##COLORS ARE JUST FOR FUN -> CHANGE AND MAKE CONSISTENT

# Define canvas and divide into two pads
c1 = ROOT.TCanvas("c", "Before and After Selection", 1200, 600)
c1.Divide(2, 1)


# Create THStack for signal and background (before selection)
stack_before = ROOT.THStack("stack_before", f"Invariant Mass of 4 leptons - Before Selection; {x_axis_titles[opt_var_index]}; Events")
sig_hist[opt_var_index].SetLineColor(ROOT.kAzure+10)
sig_hist[opt_var_index].SetMarkerColor(ROOT.kAzure+3)
sig_hist[opt_var_index].SetFillStyle(3001)
sig_hist[opt_var_index].SetLineWidth(1)

bkg_hist[opt_var_index].SetLineColor(ROOT.kMagenta+4)
bkg_hist[opt_var_index].SetFillColor(ROOT.kMagenta+4)
bkg_hist[opt_var_index].SetFillStyle(3002)
bkg_hist[opt_var_index].SetLineWidth(1)

stack_before.Add(bkg_hist[opt_var_index])
stack_before.Add(sig_hist[opt_var_index])

# Data histogram
data_hist[opt_var_index].SetLineColor(ROOT.kBlack)
data_hist[opt_var_index].SetMarkerStyle(20)
data_hist[opt_var_index].SetMarkerSize(0.8)

ymax1 = data_hist[opt_var_index].GetMaximum() * 1.3
stack_before.SetMaximum(ymax1)
# Draw on left pad
c1.cd(1)
stack_before.Draw("HIST")
data_hist[opt_var_index].Draw("E SAME")

# Legend
legend_BEF = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
legend_BEF.SetHeader(opt_var, "C")
legend_BEF.AddEntry(sig_hist[opt_var_index], "Signal")
legend_BEF.AddEntry(bkg_hist[opt_var_index], "Background")
legend_BEF.AddEntry(data_hist[opt_var_index], "Data")
legend_BEF.Draw()

# Create THStack for signal and background (after selection)
stack_after = ROOT.THStack("stack_after", f"Invariant Mass of 4 leptons - After Selection; {x_axis_titles[opt_var_index]}; Events")
sig_hist_SEL.SetLineColor(ROOT.kBlue+3)
sig_hist_SEL.SetMarkerColor(ROOT.kBlue-4)
sig_hist_SEL.SetFillStyle(3001)
sig_hist_SEL.SetLineWidth(1)

bkg_hist_SEL.SetLineColor(ROOT.kPink-7)
bkg_hist_SEL.SetFillColor(ROOT.kPink+10)
bkg_hist_SEL.SetFillStyle(3002)
bkg_hist_SEL.SetLineWidth(1)

stack_after.Add(bkg_hist_SEL)
stack_after.Add(sig_hist_SEL)

data_hist_SEL.SetLineColor(ROOT.kBlack)
data_hist_SEL.SetMarkerStyle(20)
data_hist_SEL.SetMarkerSize(0.8)

ymax2 = data_hist_SEL.GetMaximum() * 1.3
stack_after.SetMaximum(ymax2)

# Draw on right pad
c1.cd(2)
stack_after.Draw("HIST")
data_hist_SEL.Draw("E SAME")

# Legend
legend_SEL = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
legend_SEL.SetHeader(opt_var, "C")
legend_SEL.AddEntry(sig_hist_SEL, "Signal")
legend_SEL.AddEntry(bkg_hist_SEL, "Background")
legend_SEL.AddEntry(data_hist_SEL, "Data")
legend_SEL.Draw()

# Save the canvas
c1.Print(f"{opt_var}_comparison_BEFvsSEL.pdf")
c1.Clear()


## -----------------------------------------PERFORMING MAXIMUM LIKELIHOOD FITS ------------------------------------

# Initialize
mu_list = np.arange(0., 5.5, 0.01)
log_likelihood = []
graph_LL = ROOT.TGraph()

# Fill log-likelihood values and TGraph
for i, mu in enumerate(mu_list):
    ll_value = Log_Likelihood_Func(mu, sig_hist_SEL, bkg_hist_SEL, data_hist_SEL)
    log_likelihood.append(ll_value)
    graph_LL.SetPoint(i, mu, ll_value)

i_max = np.argmax(np.array(log_likelihood))
mu_hat = mu_list[i_max]

mu_intervals = compute_sigma_intervals(mu_list, log_likelihood)

print("$\hat{\mu}$ = ", mu_hat)
for sigma_label, vals in mu_intervals.items():
    print(f"{sigma_label} interval: [{vals['lower']:.3f}, {vals['upper']:.3f}], "
          f"uncertainty ±{vals['half_width']:.3f}")

# Optional: Draw the graph
c3 = ROOT.TCanvas("c", "Log Likelihood vs #mu", 800, 600)
graph_LL.SetTitle("Log-Likelihood vs. Signal Strength #mu")
graph_LL.GetXaxis().SetTitle("#mu")
graph_LL.GetYaxis().SetTitle("Log-Likelihood")
graph_LL.SetLineColor(ROOT.kBlue + 1)
graph_LL.SetLineWidth(2)
graph_LL.Draw("AL")  # 'A' to draw axes, 'L' for line

# Draw vertical line at mu_hat
y_min = graph_LL.GetYaxis().GetXmin()
y_max = graph_LL.GetYaxis().GetXmax()
line = ROOT.TLine(mu_hat, y_min, mu_hat, y_max)
line.SetLineColor(ROOT.kMagenta)
line.SetLineStyle(2)  # dashed line
line.SetLineWidth(2)
line.Draw("same")

# Optional: Add legend or label for mu_hat
label = ROOT.TLatex(mu_hat, y_max - 0.2 * (y_max - y_min), f"#hat{{#mu}} = {mu_hat:.2f}")
label.SetTextColor(ROOT.kMagenta+3)
label.SetTextSize(0.04)
label.SetTextAlign(21)
label.Draw()

# Define color for 1σ region
color_1sigma = ROOT.kGreen + 1

# Extract 1σ interval
vals = mu_intervals["1σ"]
lower = vals["lower"]
upper = vals["upper"]

# Safety check
if lower is not None and upper is not None:
    # Draw vertical lines
    line_lower = ROOT.TLine(lower, y_min, lower, y_max)
    line_lower.SetLineColor(color_1sigma)
    line_lower.SetLineStyle(2)
    line_lower.SetLineWidth(1)
    line_lower.Draw("same")

    line_upper = ROOT.TLine(upper, y_min, upper, y_max)
    line_upper.SetLineColor(color_1sigma)
    line_upper.SetLineStyle(2)
    line_upper.SetLineWidth(1)
    line_upper.Draw("same")

    # Draw shaded region
    box = ROOT.TBox(lower, y_min, upper, y_max)
    box.SetFillColorAlpha(color_1sigma, 0.15)
    box.Draw("same")

    # Add label
    label_text = "68% CL (1#sigma interval)"
    label_x = upper * 1.1
    label_y = y_max - 0.4 * (y_max - y_min)

    text = ROOT.TLatex(label_x, label_y, label_text)
    text.SetTextAlign(21)  # center alignment
    text.SetTextColor(color_1sigma + 3)
    text.SetTextSize(0.035)
    text.Draw("same")


c3.Update()

c3.Print("LogLikelihood_vs_mu.pdf")

##For somo reason, writing #mu yields a proprotional sign?? LOOK INTO THIS!!
# Happens only when printing into pdf, but, otherwise, looses lots of quality...





















