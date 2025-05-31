import ROOT
import os

##This code converts the needed ROOT files onto .csv files. 
# Such is done because we want to optimize our signal selection using Machine Learning, hence, we must have our data in .csv files.

##Function to ransform each individual histogram onto a CSV file
def histogram_to_csv(histogram, csv_filename):

    with open(csv_filename, "w") as csv_file:
        csv_file.write("BinIndex,BinContent\n")

        for i in range(1, histogram.GetNbinsX() + 1):
            bin_index = i
            bin_content = histogram.GetBinContent(i)
            csv_file.write(f"{bin_index},{bin_content}\n")

    print(f"Saved: {csv_filename}")

def Convert_Root_Files_Func(input_dir, output_dir, histogram_names):

    os.makedirs(output_dir, exist_ok=True)

    for root_file_name in os.listdir(input_dir):

        if not root_file_name.endswith(".root"):
            continue

        root_path = os.path.join(input_dir, root_file_name)
        root_file = ROOT.TFile(root_path)

        if not root_file or root_file.IsZombie():
            print(f"Error reading {root_path}")
            continue

        # Strip .root extension for naming
        file_prefix = os.path.splitext(root_file_name)[0]

        for histogram_name in histogram_names:

            histogram = root_file.Get(histogram_name)

            if histogram:
                csv_name = f"{file_prefix}_{histogram_name}.csv"
                csv_path = os.path.join(output_dir, csv_name)
                histogram_to_csv(histogram, csv_path)

            else:
                print(f"Warning: '{histogram_name}' not found in {root_file_name}")

        root_file.Close()

if __name__ == "__main__":
    
    input_directory = "./Data_MC_Samples"       # Folder containing all .root files
    output_directory = "./csv_output"    # Folder to save .csv files

    # List of histograms to extract from each file

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

    Convert_Root_Files_Func(input_directory, output_directory, histogram_names)
