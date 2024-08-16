# Bacterial_image_analysis
## Macros and methods to analyse bacterial bioimages

This repository describes the macros used in the manuscript **Transertion and cell geometry organize the Escherichia coli nucleoid during rapid growth** found [here](https://doi.org/10.1101/2023.10.16.562172).

## Python code to analyse the distance between the bacterial membrane and DNA

1. Clone the repository:
Open up the terminal and paste the following commands
````
git clone https://github.com/CKSpahn/Bacterial_image_analysis.git
cd Bacterial_image_analysis
mamba env create -f environment.yml
mamba activate bacteria_ia
````

2. Open the file `dna-membrane-dist.py` inside the Python folder and change the following paths accordingly:

````
root_folder = '/Users/user1/Documents/data/DNA_membrane_distance_diff_treatments_2024-08-05'
results_folder = '/Users/user1/Documents/results/'
````


