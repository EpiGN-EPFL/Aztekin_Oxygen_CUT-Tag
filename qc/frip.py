import numpy as np
import pandas as pd
import os
import glob
import argparse

def read_dfs_in_directory(directory, pattern="summary", verbose=False):
    # Check if the directory exists
    if not os.path.isdir(directory):
        print(f"The directory {directory} does not exist.")
        return
    
    dfs=[]
    # Read all files in the directory as dataframes 
    for filename in os.listdir(directory):
        if pattern in filename:
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                if verbose:
                    print(f"Reading file: {file_path}")
                with open(file_path, 'r') as file:
                    dfs.append(pd.read_csv(file, sep='\t', index_col=0))
    
    # Concat to the final df
    df=pd.concat(dfs, axis=1)

    return(df)

def calculate_FRiP(df):
    # Calculate total reads
    df.loc['Total_Read']=df.sum(axis=0)

    # Calculate FRiP (Fraction of Reads in Peaks)
    df.loc['FRiP']=df.loc['Assigned']/df.loc['Total_Read']

    return(df)
    
def parse_args():
    # Create the parser
    parser = argparse.ArgumentParser(description="Read featureCounts from a specified directory and calculate the FRiP score.")

    # Add the argument
    parser.add_argument('-i', '--input_dir', type=str, required=True, 
                        help="The path to the directory containing files to read.")
    parser.add_argument('-p', '--pattern', type=str, required=False, 
                        default="summary", help="Patten of the files to read in input")
    parser.add_argument('-v', '--verbose', type=bool, required=False, default=False)
    parser.add_argument('-o', '--output_file', type=str, required=True)

    # Parse the arguments
    args = parser.parse_args()

    return(args)

def main():
    # Read summary from featureCounts
    args=parse_args()
    print(f"Reading files from {args.input_dir}")
    df = read_dfs_in_directory(args.input_dir, args.pattern)
    df = calculate_FRiP(df)
    df.to_csv(args.output_file, sep='\t')
    print(f"Saving result to {args.output_file}")

if __name__ == "__main__":
    main()