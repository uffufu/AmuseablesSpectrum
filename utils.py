import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# ===== cube reading and spectrum extraction =====
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
# ===== file handling =====
import os
from glob import glob
import tkinter as tk
from tkinter import filedialog
from datetime import datetime

def select_file():
    # ===== user selects a single FITS file. =====
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(title="Select FITS file", filetypes=[("FITS files", "*.fits")])

    return file_path

def select_folder_and_find_fits():
    # ===== user selects a folder, and the script finds all FITS files in that folder. =====
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(title="Select folder")
    if not folder_path:
        print("No folder selected.")
        return [], None

    fits_files = glob(os.path.join(folder_path, "*.fits"))
    if not fits_files:
        print("No FITS files found in the selected folder.")
        return [], None

    folder_name = os.path.basename(folder_path)
    print("FITS files found：")
    for i, f in enumerate(fits_files):
        print(f"{i}: {f}")

    return fits_files, folder_name

def read_fits_file(file_name):
    # ===== read the FITS file and return the data and header. =====
    cube_data = fits.getdata(file_name)
    header = fits.getheader(file_name)
    n_channels = cube_data.shape[0]
    print(f"total channels: {n_channels}")
    base = os.path.splitext(os.path.basename(file_name))[0]

    return cube_data, header, n_channels, base

def get_spectrum(cube_data):
    # ===== calculate the spectrum from the cube data. =====
    spectrum = np.nanmean(cube_data, axis=(1, 2)) 

    return spectrum

def exclude_edge_channels(n_channels, exclude_fraction=1/8):
    # ===== exclude the edge channels based on the specified fraction. =====
    start_idx = int(n_channels * exclude_fraction)
    end_idx = int(n_channels * (1 - exclude_fraction))
    print(f"channels used: {start_idx+1} to {end_idx}") # +1 because the end index start from 0
    
    return start_idx, end_idx

def filter_out_channels(spectrum, start_idx, end_idx):
    # ===== roughly filter out channels by checking if the signal is above RMS of spectrum intensity =====

    # calculate the RMS of the original spectrum
    original_rms = np.sqrt(np.nanmean(spectrum[start_idx:end_idx]**2))
    # check if the signal above the RMS is continuous for at least 3 channels
    valid_channels = []
    for i in range(start_idx, end_idx - 2):# -2 is because we check 3 channels
        if all(spectrum[i:i+3] > original_rms):
            valid_channels.extend([i, i+1, i+2])
        
    valid_channels = sorted(set(valid_channels))
    ranges = []
    if valid_channels:
        start = valid_channels[0]
        for i in range(1, len(valid_channels)):
            if valid_channels[i] != valid_channels[i - 1] + 1:  # if not continuous
                ranges.append((start, valid_channels[i - 1]))  # add range
                start = valid_channels[i]  # update start point
        ranges.append((start, valid_channels[-1]))  # add the last range
    print(f"number of channels with signals: {len(valid_channels)}, channels with signals: {ranges}")
    # plot
    plt.figure(figsize=(15, 3))
    plt.plot(spectrum, label='Spectrum')
    plt.axhline(y=original_rms, color='green', linestyle='--', label=f'Rough Threshold = {original_rms:.4f}')
    for ch in valid_channels:
        plt.axvline(x=ch, color='yellow', alpha=0.5, label='Filtered Channels' if ch == valid_channels[0] else None)# just for label
    
    plt.title("Original Spectrum with Threshold")
    plt.xlabel("Channel")
    plt.ylabel("Intensity")
    plt.legend()
    plt.show()

    return valid_channels, ranges

def masking(cube_data, valid_channels, start_idx, end_idx):
    # ===== create masked cube data =====

    # calculate RMS of each noise channel and take the average
    noise_channels = [i for i in range(start_idx, end_idx+1) if i not in valid_channels]
    noise_rms = [np.sqrt(np.nanmean(cube_data[ch, :, :]**2)) for ch in noise_channels]
    rms = np.mean(noise_rms)

    # only process valid channels
    valid_cube = cube_data[valid_channels, :, :]

    # apply 3𝜎 mask
    masked_valid_cube = np.where(valid_cube >= 3 * rms, valid_cube, np.nan)
    
    return masked_valid_cube, rms

def moment_maps(masked_valid_cube, header, base):
    # ===== use SpectralCube to calculate Moment Maps =====

    wcs = WCS(header)
    cube = SpectralCube(data=masked_valid_cube, wcs=wcs)
    
    moment0 = cube.moment(order=0)
    moment1 = cube.moment(order=1)
    moment2 = cube.moment(order=2)
    
    # convert units to be consistent with CARTA units
    if moment0.unit.is_equivalent(u.Jy / u.beam * u.m / u.s):
        moment0 = moment0.to(u.Jy / u.beam * u.km / u.s)
    if moment1.unit.is_equivalent(u.m / u.s):
        moment1 = moment1.to(u.km / u.s)
    if moment2.unit.is_equivalent((u.km / u.s)**2):
        moment2_disp = (moment2 ** 0.5).to(u.km / u.s) # dispersion

    # Generate a timestamp
    timestamp = datetime.now().strftime("%m%d_%H%M%S")

    # Save Moment Maps as FITS files
    moment0.write(f"{base}_moment0_{timestamp}.fits", overwrite=True)
    moment1.write(f"{base}_moment1_{timestamp}.fits", overwrite=True)
    moment2_disp.write(f"{base}_moment2_{timestamp}.fits", overwrite=True)
    print(f"Moment Maps saved as{base}_moment0/1/2_{timestamp}.fits")

def write_to_excel(file_name, rms, ranges, excel_path, folder_name=False):
    # ===== write the results to an Excel file, sheet name based on file or folder name =====
    
    output_data = {
        "filename": [file_name],
        "time": [datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        "rms": [rms],
    }
    # write ranges to output_data
    for idx, (start, end) in enumerate(ranges):
        output_data[f"range_{idx + 1}_start"] = [start]
        output_data[f"range_{idx + 1}_end"] = [end]

    df = pd.DataFrame(output_data)

    # worksheet name based on file or folder name
    if folder_name and os.path.dirname(file_name):
        sheet_name = os.path.basename(os.path.dirname(file_name))
    else:
        sheet_name = os.path.splitext(os.path.basename(file_name))[0]
    # limit sheet name length to 31 characters
    sheet_name = sheet_name[:31]

    # check if the directory exists, create it if not
    
    if os.path.dirname(excel_path):
        os.makedirs(os.path.dirname(excel_path), exist_ok=True)

    # if Excel file already exists, read existing data and append
    if os.path.exists(excel_path):
        with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a" ,if_sheet_exists="overlay") as writer:
            try:
                existing_df = pd.read_excel(excel_path, sheet_name=sheet_name)
                df = pd.concat([existing_df, df], ignore_index=True)
            except ValueError:
                # if sheet does not exist, write directly
                pass
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        # if the file does not exist, create a new one
        with pd.ExcelWriter(excel_path, engine="openpyxl", mode="w") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Write in Excel: {excel_path}，Sheet: {sheet_name}")