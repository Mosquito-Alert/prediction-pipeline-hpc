import argparse
import rasterio
import numpy as np

def classify_selected(value_array, nodata_value):
    """
    Classify raster values into categories with numeric labels:
        - 1: discont_urban_fabric (value == 2)
        - 2: other_artificial (values 3, 5–9)
        - 3: agricultural (values 12–22)
        - 0: other (all other values)
    """

    label_map = {
        'other': 0,
        'discont_urban_fabric': 1,
        'other_artificial': 2,
        'agricultural': 3,
    }

    # Define output nodata value (must not conflict with 0–3)
    output_nodata = -128

    if nodata_value:
        nodata_mask = value_array == nodata_value
    else:
        nodata_mask = np.zeros(value_array.shape, dtype=bool)

    # Initialize classification
    classified = np.zeros(value_array.shape, dtype=np.int8)  # Default 0: other

    # Apply classification only to valid (non-nodata) pixels
    valid_data = ~nodata_mask
    classified[(value_array == 2) & valid_data] = label_map['discont_urban_fabric']
    classified[(np.isin(value_array, [3, 5, 6, 7, 8, 9])) & valid_data] = label_map['other_artificial']
    classified[(np.isin(value_array, range(12, 23))) & valid_data] = label_map['agricultural']

    # Apply nodata to output
    classified[nodata_mask] = output_nodata

    return classified, output_nodata

def process_tif(input_file, output_file):
    with rasterio.open(input_file) as src:
        data = src.read(1)  # Read first band
        profile = src.profile

        classified_array, output_nodata = classify_selected(data, nodata_value=src.nodata)

        # Update profile for single-band int8 output
        profile.update(dtype=rasterio.int8, count=1, nodata=output_nodata, compress='deflate')  # compress optional but good for size

        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(classified_array, 1)

def main():
    parser = argparse.ArgumentParser(description="Classify selected CLC categories from a GeoTIFF file.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to input .tif file")
    parser.add_argument("-o", "--output_file", required=True, help="Path to output .tif file")
    args = parser.parse_args()

    process_tif(input_file=args.input_file, output_file=args.output_file)

if __name__ == "__main__":
    main()