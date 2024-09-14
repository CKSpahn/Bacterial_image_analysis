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

## Fiji macros

### Overview

The following table specifies the individual macros used in the manuscript

| Macro | Purpose | Input | Input format | #Channels | Output | Comments | Link |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M1 | Extracting of the experimental PSF | z-stack | TIFF | 1 | Individual PSFs, Average PSF | - | [Macro M1](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M1%20-%20Average%20PSF%20from%20beadstack.ijm) |
| M2 | Smoothing of cell outlines for selection of individual bacteria | Single- or multi-channel SMLM images | TIFF | 1+ (in this study: 3) | Smoothed membrane image, opened in Fiji | - | [Macro M2](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M2%20-%20Smoothing%20of%20cell%20outlines%20in%20NR%20PAINT%20images.ijm) |
| M3 | Rotation, straightening of bacterial cells | CLSM image, ROIs | TIFF | 3 | Single bacteria as individual TIFF images | ROIs need to be opened in Fiji, adjust scaling if required | [Macro M3](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M3%20-%20Rotation%2C%20alignment%20and%20cell%20straightening%20-%20CLSM.ijm) |
| M4 | Normalization of straightened bacteria | Output of Macro M3 | TIFF | 3 | Normalized cells as individual TIFF image | Membrane: Channel 1 | [Macro M4](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M4%20-%20Cell%20normalization%20-%20CLSM.ijm) |
| M5 | Determination of relative nucleoid expansion | Output  Macro 3 | TIFF | 3 | .txt file summaring cell and nucleoid length, as well as the relative nucleoid expansion | Membrane: Channel 1, DNA: Channel 2 | [Macro M5](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M5%20-%20%20Determination%20of%20the%20relative%20nucleoid%20length%20expansion%20(RNLE)%20-%20CLSM.ijm) |
| M6 | Determination of the relative MreB distribution | Population average image | TIFF | 3 | .txt file with MreB intensities at poles and cell cylinder, ROIs | Membrane: Channel 1, DNA: Channel 2, MreB: Channel 3  | [Macro M6](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M6%20-%20Determination%20of%20the%20relative%20MreB%20distribution%20-%20CLSM.ijm) |
| M7 | Erosion analysis for radial intensity distribution analysis | SMLM image | TIFF | 3 | .txt files with intensity distributions for each cell, Eroded ROIs as .zip, overlay PNG | Channel order can be defined | [Macro M7](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M7%20-%20Determination%20of%20radial%20intensity%20distributions%20using%20erosion%20analysis.ijm) |
| M8 | Plotting intensities along cell perimeter | SMLM images, ROIs | TIFF | 2+ | Intensity traces of DNA and MreB channels | Image and ROIs need to be opened in Fiji | [Macro M8](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M8%20-%20Plot%20profiles%20along%20perimenter.ijm) |
| M9 | Simulation of 3-color images for circular x-corr | - | - | 3 | Simulated 3-color image and ROIs | Simulate points within a circular cell with varying number, shift and radial position | [Macro M9](https://github.com/CKSpahn/Bacterial_image_analysis/blob/main/Macros/Macro_M9%20-%20Simulate%20images%20for%20circular%20cross-correlation.ijm) |
