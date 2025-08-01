# AmuseablesSpectrum

Amuseables Spectrum is tuned into cosmic waves from across the universe and beyond!
From distant baby sun to dancing disk, the antenna arrays are always listening, picking up whispers from momentum. 
Now looking for a co-pilot for a new wave of celestial adventures. Ready to cruise the cosmic frequencies and see what the universe is broadcasting today?
Pack your data cube and let's ride the radio waves to wherever the stars want to take us!

Automated spectral data processing tool that supports single or batch FITS cubes, generates moment map FITS files, and outputs analysis results to Excel.

Let's think about the cost of building a space telescope.
YOUR DATA CUBES ARE SO EXPENSIVE🤑!
PLEASE CHECK THE ORIGINAL DATACUBES WITH YOUR EYES JUST TO BE SAFE!

## Features 
- Support single FITS cube or batch processing of all FITS cubes in a folder
- Detect frequency range automatically and estimate RMS
- Generate moment 0/1/2 FITS file
- RMS and channel used written into Excel（one sheet for each file）

## Working Principles

The program automatically extracts mean spectrum from FITS cubes, identifies valid channels with signal above noise level, applies 3σ threshold masking, and generates moment maps (integrated intensity, velocity field, and velocity dispersion) for astronomical analysis.

## Installation

1. Python 3.9（or compatible version）recommended
2. install required packages：
   ```bash
   pip install -r requirements.txt
   ```

## Instructions

### 1. Run main.py directly
```bash
python main.py
```
Or Directly Run All Cells of Jupyter Notebook `AmuseablesSpectrum.ipynb` . 

### 2. Mode selection
- enter `f`：for select single FITS cube manually
- enter `d`：select a folder, automatically process all FITS cubes in the folder

### 3. Output
- moment map is saved to the execution directory (or specified folder)
- Statistical data is saved to `moment_info.xlsx`

## Structure

```
AmuseablesSpectrum/
├── AmuseablesSpectrum.ipynb    # Include main program and tool functions in one Jupyter Notebook
├── main.py           # Main program
├── utils.py          # Tool function
├── requirements.txt  # Packages list
├── README.md         # Instruction manual
```

## Acknowledgement

Part of the programming and documentation for this project was completed with the help of AI tools.