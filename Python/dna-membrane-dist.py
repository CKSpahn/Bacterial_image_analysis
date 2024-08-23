import numpy as np
import os
import pandas as pd
import re

import matplotlib.pyplot as plt
import seaborn as sns

root_folder = '/Users/esti/Documents/PROYECTOS/NUCLEOID-ECOLI/data/DNA_membrane_distance_diff_treatments_2024-08-05'
results_folder = '/Users/esti/Documents/PROYECTOS/NUCLEOID-ECOLI/'



def extract_data(file_path, intensity_type):
    """
    Extracts the distance between the maximum intensities of Membrane and DNA from a CSV file.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        float: The distance between the maximum intensities of Membrane and DNA.
    """
    # Read the text file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Find the row with maximum intensity for Membrane_intensity and DNA_intensity
    membrane_max = df['Membrane_intensity'].idxmax()
    if intensity_type == 'normalized' or intensity_type == 'normalised':
        mask = df['DNA_intensity'] > (1/np.exp(1))
    else:
        th = np.max(df['DNA_intensity'])/np.exp(1)
        mask = df['DNA_intensity'] > th
        print(th)
    dna_max = mask.idxmax()
    # Calculate the distance between Membrane and DNA max intensities
    distance = abs(df.loc[membrane_max, 'distance'] - df.loc[dna_max, 'distance'])

    return distance


def extract_information_from_filename(filename):
    """
    Extracts cell and ROI information from a filename.

    Args:
        filename (str): The name of the file.

    Returns:
        Tuple[Optional[str], Optional[str], str]: A tuple containing the cell and ROI information.
            The cell information is the part of the filename that starts with 'm', or None if not found.
            The ROI information is the part of the filename that starts with 'ROI', or None if not found.
            The intensity_type is the last part of the filename that contains the word 'normalized', or 'raw' if not found.
    """

    # Extract cell information
    cell_match = re.search(r'm(\d+)', filename)
    cell = f"m{cell_match.group(1)}" if cell_match else None

    # Extract ROI information
    roi_match = re.search(r'ROI_(\d+)', filename)
    roi = f"ROI_{roi_match.group(1)}" if roi_match else None

    # Determine the intensity type based on the last part of the filename
    # Split the filename into parts
    name_parts = filename.split('_')
    intensity_type = "normalized" if "normalized" in name_parts[-1] else "raw"

    return cell, roi, intensity_type


def process_folders(root_folder):
    """
    Process all txt files in the given root folder and its subdirectories.
    Returns a pandas DataFrame containing information about the distances between the membrane and DNA maximum intensities.
    Args:
        root_folder (str): The path to the root folder to search for txt files.
    Returns:
        pandas.DataFrame: A DataFrame containing the following columns:
            - condition (str): The name of the condition folder.
            - timepoint (str): The name of the timepoint folder.
            - file (str): The name of the txt file.
            - distance (float): The extracted distance from the txt file.
            - cell (str): The extracted cell information from the filename.
            - roi (str): The extracted ROI information from the filename.
            - intensity_type (str): Whether the analysis comes from normalised intensities or raw intensity values.
    """
    results = []

    # Walk through all subdirectories
    for condition in os.listdir(root_folder):
        condition_path = os.path.join(root_folder, condition)
        if os.path.isdir(condition_path):
            for timepoint in os.listdir(condition_path):
                timepoint_path = os.path.join(condition_path, timepoint)
                if os.path.isdir(timepoint_path):
                    # Using os.listdir to find all txt files
                    for file in os.listdir(timepoint_path):
                        if file.endswith('.txt'):
                            file_path = os.path.join(timepoint_path, file)
                            # Exclude .DScore files and hidden files
                            if not file.endswith('.DScore') and not file.startswith('.'):
                                try:
                                    cell, roi, intensity_type = extract_information_from_filename(file)
                                    distance = extract_data(file_path, intensity_type)
                                    results.append({
                                        'condition': condition,
                                        'timepoint': timepoint,
                                        'file': os.path.basename(file),
                                        'distance': distance,
                                        'cell': cell,
                                        'roi': roi,
                                        'type': intensity_type
                                    })
                                except Exception as e:
                                    print(f"Error processing file {file}: {str(e)}")

    return pd.DataFrame(results)


# Extract distances from all the conditions
results_df = process_folders(root_folder)

# Save results to a CSV file
results_df.to_csv(os.path.join(results_folder, 'results.csv'), index=False)
print("Data extraction complete. Results saved to results.csv")

### STATISTICS AND PLOTTING ###

# Read the CSV file
df = pd.read_csv(os.path.join(results_folder, 'results.csv'))

# Convert timepoint to numeric, assuming format like '10_min'
df['timepoint'] = df['timepoint'].str.extract('(\d+)').astype(int)

# Separate normalized and raw data
df_normalized = df[df['type'] == 'normalized']
df_raw = df[df['type'] == 'raw']

# Calculate average distance for each cell at each time point for normalized and raw data
avg_distances_normalized = df_normalized.groupby(['condition', 'timepoint', 'cell'])['distance'].mean().reset_index()
avg_distances_raw = df_raw.groupby(['condition', 'timepoint', 'cell'])['distance'].mean().reset_index()

# Calculate the mean and standard deviation for each condition and time point
population_stats_normalized = avg_distances_normalized.groupby(['condition', 'timepoint']).agg({
    'distance': ['mean', 'std']
}).reset_index()
population_stats_normalized.columns = ['condition', 'timepoint', 'mean_distance', 'std_distance']

population_stats_raw = avg_distances_raw.groupby(['condition', 'timepoint']).agg({
    'distance': ['mean', 'std']
}).reset_index()
population_stats_raw.columns = ['condition', 'timepoint', 'mean_distance', 'std_distance']

# Sort the dataframes by condition and timepoint
population_stats_normalized = population_stats_normalized.sort_values(['condition', 'timepoint'])
population_stats_raw = population_stats_raw.sort_values(['condition', 'timepoint'])

# Create the plot
fig = plt.figure(figsize=(18, 7))
sns.set_style("whitegrid")
plt.subplot( 1, 2, 1)
# Plot each condition for normalized data
for condition in population_stats_normalized['condition'].unique():
    condition_data = population_stats_normalized[population_stats_normalized['condition'] == condition]
    plt.errorbar(condition_data['timepoint'], condition_data['mean_distance'],
                 yerr=condition_data['std_distance'], capsize=3,
                 marker='o', linestyle='-')
plt.xlabel('Time Point (min)')
plt.ylabel('Average Distance')
plt.title('Temporal Curves of Average Distances by Condition (Normalized)')
#plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.subplot(1, 2, 2)
# Plot each condition for raw data
for condition in population_stats_raw['condition'].unique():
    condition_data = population_stats_raw[population_stats_raw['condition'] == condition]
    plt.errorbar(condition_data['timepoint'], condition_data['mean_distance'],
                 yerr=condition_data['std_distance'], capsize=3,
                 label=f'{condition} (Raw)', marker='o', linestyle='-')
plt.xlabel('Time Point (min)')
plt.ylabel('Average Distance')
plt.title('Temporal Curves of Average Distances by Condition (Raw)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

fig.tight_layout()
fig.savefig(os.path.join(results_folder, 'temporal_curves.png'))
plt.show()

# Print summary statistics and save data
print("Normalized Data Statistics:")
population_stats_normalized.to_csv(os.path.join(results_folder, 'stats_normalised_results.csv'), index=False)
print(population_stats_normalized)
print("\nRaw Data Statistics:")
population_stats_raw.to_csv(os.path.join(results_folder, 'stats_raw_results.csv'), index=False)
print(population_stats_raw)