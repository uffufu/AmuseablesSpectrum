from utils import (
    select_file, 
    select_folder_and_find_fits, 
    read_fits_file, get_spectrum, 
    exclude_edge_channels, 
    filter_out_channels, 
    masking, moment_maps, 
    write_to_excel)
'''main workflow'''
fraction = 1/8  # default fraction of channels to exclude from edges
excel_path = "moment_info.xlsx"  # default path for Excel file

mode = input("Enter ‘f’ to select manually, or ‘d’ to select folder and automatically search for FITS files:")
if mode == 'f':
    file_list = [select_file()]
    folder_name = None
elif mode == 'd':
    file_list, folder_name = select_folder_and_find_fits()
else:
    raise ValueError("Please enter ‘f’ or 'd'.")

print("Files used:")
for file_name in file_list:
    print(file_name)

    cube_data, header, n_channels, base = read_fits_file(file_name)
    spectrum = get_spectrum(cube_data)

    start_idx, end_idx = exclude_edge_channels(n_channels, fraction)
    valid_channels, ranges = filter_out_channels(spectrum, start_idx, end_idx)

    masked_valid_cube, rms = masking(cube_data, valid_channels, start_idx, end_idx)
    
    moment_maps(masked_valid_cube, header, base)

    write_to_excel(base, rms, ranges, excel_path)