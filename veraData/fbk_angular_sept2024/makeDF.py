import os
import pandas as pd

# Directory containing the CSV files
directory = '/Users/austindestcroix/Desktop/school_holders/projects/vera/daaps/veraData/fbk_angular_sept2024/'

# List to hold individual DataFrames
dfs = []

# Loop through all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path)
        
        # find the number before 'nm' in filepath
        nm_index = filename.find('nm')
        wl = float(filename[nm_index-3:nm_index])
        print(wl)
        
        # df['source_file'] = filename  # Add a column to identify the source file
        df['wavelength (nm)'] = wl
        
        temp_index = filename.find('C')
        temp = 273+float(filename[nm_index+3:temp_index])
        print(temp)
        df['temperature (K)'] = temp
        
        df['voltage (V)'] = 6.0
        
        dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_df = pd.concat(dfs, ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file = os.path.join(directory, 'combined_long_format.csv')
combined_df.to_csv(output_file, index=False)

print(f"Combined CSV saved to {output_file}")