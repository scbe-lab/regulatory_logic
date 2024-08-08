import argparse
import pandas as pd
import numpy as np

def get_gie_stain(gc_content):
    quantiles = [0, 0.25, 0.5, 0.75, 0.9]
    q_values = np.quantile(gc_content, quantiles)
    if gc_content >= q_values[4]:
        return "gpos100"
    elif gc_content >= q_values[3]:
        return "gpos75"
    elif gc_content >= q_values[2]:
        return "gpos50"
    elif gc_content >= q_values[1]:
        return "gpos25"
    else:
        return "gneg"

def main(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", comment="#", header=None, names=["chrom", "start", "end", "name", "pct_at", "pct_gc", "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len"])
    df["gieStain"] = df["pct_at"].apply(get_gie_stain)

    with open(output_file, "w") as outfile:
        for _, row in df.iterrows():
            outfile.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t{row['gieStain']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GC content to cytoband format.")
    parser.add_argument("input_file", help="Path to the input gc_content.txt file")
    parser.add_argument("output_file", help="Path to the output cytoband.txt file")
    args = parser.parse_args()

    main(args.input_file, args.output_file)
