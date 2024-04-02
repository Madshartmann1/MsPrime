import sys
import pandas as pd
import gzip
from tqdm import tqdm
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

def print_with_timestamp(message):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"{timestamp}: {message}")


def read_vcf(vcf_file):
    print_with_timestamp("Reading the VCF file")
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    with open_func(vcf_file, 'rt') as f:
        header_line_num = None
        for line_num, line in enumerate(f):
            if line.startswith('#CHROM'):
                header_line_num = line_num
                break
        if header_line_num is not None:
            f.seek(0)
            vcf_df = pd.read_csv(f, sep='\t', header=header_line_num)
    print_with_timestamp("Finished reading the VCF file")
    return vcf_df


def read_intervals(intervals_file):
    """
    Parses the intervals file, organizing data by individual and chromosome type, including
    the intervals that are on the same line as the M or P chromosome designation.
    
    Args:
        intervals_file (str): Path to the intervals file.
    
    Returns:
        list of tuples: List containing tuples of (individual, chrom, start_pos, end_pos, barcode).
    """

    print_with_timestamp(f"Reading intervals from {intervals_file}...")
    
    intervals = []
    current_individual = None
    current_chrom = None


    # Adjust to handle gzipped intervals file
    open_func = gzip.open if intervals_file.endswith('.gz') else open
    
    with open_func(intervals_file, 'rt') as file:  # 'rt' mode for reading as text
        line = None  # Initialize line before the loop
        for line in file:
            line = line.strip()
            # Check for individual designation
            if line.startswith("Individual"):
                current_individual = line.split()[1]
            else:
                # Split the line, accounting for possible M or P designations with data
                parts = line.split()
                if parts and parts[0] in ("M", "P"):
                    current_chrom = parts[0]
                    # Check if there's interval data following the chromosome designation on the same line
                    if len(parts) > 1:
                        # Combine the parts back and split by comma for interval data
                        interval_data = ' '.join(parts[1:]).split(', ')
                        if len(interval_data) == 3:
                            start_pos, end_pos, barcode = interval_data
                            intervals.append((current_individual, current_chrom, int(start_pos), int(end_pos), barcode))
                elif parts:  # Processing interval lines without M or P designation
                    interval_data = line.split(', ')
                    if len(interval_data) == 3:
                        start_pos, end_pos, barcode = interval_data
                        intervals.append((current_individual, current_chrom, int(start_pos), int(end_pos), barcode))
    
    print_with_timestamp("Finished reading intervals.")
    return intervals


def map_barcode_to_vcf_column(barcode):
    """
    Maps a barcode to the corresponding VCF column index.

    Args:
        barcode (str): The barcode string (e.g., 'KAR#57_P').

    Returns:
        int: The column index in the VCF DataFrame corresponding to the barcode.
    """

    # Split the barcode into its components
    parts = barcode.split('#')
    prefix = parts[0]
    num_suffix = parts[1].split('_')
    num = int(num_suffix[0])

    # Define the base indices mapping
    base_indices = {"CEU": 0, "YRI": 500, "CHB": 1000, "KAR": 1500}

    # Calculate the tsk_# index
    tsk_index = base_indices[prefix] + num

    # Assuming the VCF DataFrame has a direct mapping of tsk_index to column names,
    # like 'tsk_0', 'tsk_1', ..., adjust as necessary for your specific VCF format.
    # The +9 offset is based on VCF format specifics, adjust as needed.
    vcf_col_name = f"tsk_{tsk_index}"

    return vcf_col_name




def process_vcf_intervals(vcf_df, intervals):
    """
    Processes the VCF DataFrame based on intervals, extracting and concatenating allele information
    for M and P chromosomes for each individual.

    Args:
        vcf_df (pd.DataFrame): DataFrame containing the VCF data.
        intervals (list): A list of tuples with interval data (individual, chrom, start_pos, end_pos, barcode).

    Returns:
        dict: A dictionary where keys are individual IDs and values are strings of concatenated alleles in "M|P" format.
    """
    # Initialize dictionaries to temporarily hold M and P allele information
    alleles_temp = {'M': {}, 'P': {}}
    
    # Initialize the final dictionary for processed alleles
    processed_alleles = {}

    
    for individual, chrom, start_pos, end_pos, barcode in tqdm(intervals, desc="Processing intervals"):
       # Filter the DataFrame for rows where POS is within the interval
        filtered_df = vcf_df[(vcf_df['POS'] >= start_pos) & (vcf_df['POS'] <= end_pos)]
        
        # Use the mapping function to get the corresponding VCF column name
        vcf_col_name = map_barcode_to_vcf_column(barcode)


        # Initialize a list to hold alleles for this interval
        alleles_list = []

        # Iterate through the filtered DataFrame and extract alleles
        for _, row in filtered_df.iterrows():
            alleles = row[vcf_col_name].split('|')
            
            # Determine lineage from barcode ('M' for maternal, 'P' for paternal)
            lineage = 'M' if 'M' in barcode else 'P'
            
            # Choose allele based on lineage, not chrom
            allele = alleles[0] if lineage == 'M' else alleles[1]
            alleles_list.append(allele)
    
              
            # Accumulate alleles for both M and P data
        if individual not in alleles_temp[chrom]:
            alleles_temp[chrom][individual] = alleles_list
        else:
            # Extend the existing allele list for this individual and chromosome
            alleles_temp[chrom][individual].extend(alleles_list)
    
    # After processing all intervals, concatenate M and P data for each individual
    for individual in set(ind for ind, _, _, _, _ in intervals):
        m_alleles = alleles_temp['M'].get(individual, [])
        p_alleles = alleles_temp['P'].get(individual, [])
        # Concatenate M and P alleles
        concatenated_alleles = [f"{m}|{p}" for m, p in zip(m_alleles, p_alleles)]
        processed_alleles[individual] = concatenated_alleles

    return processed_alleles




def save_vcf_with_headers(vcf_df, output_file):
    open_func = gzip.open if output_file.endswith('.gz') else open
    with open_func(output_file, 'wt') as f:
        # Write the VCF standard header lines
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=tskit 0.5.4\n")
        f.write("##FILTER=<ID=PASS,Description=All filters passed>\n")
        f.write("##contig=<ID=1,length=248387328>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
        
        # Fixed VCF headers for the first 9 columns
        fixed_headers = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        
        # Extracting additional headers (individual identifiers) and sorting them
        additional_headers = [col for col in vcf_df.columns if col not in fixed_headers]
        sorted_additional_headers = sorted(additional_headers, key=lambda x: int(x))  # Assuming identifiers are numeric
        
        # Combine fixed headers with sorted additional headers
        all_headers = fixed_headers + sorted_additional_headers
        
        # Write the full header line to the file
        f.write('\t'.join(all_headers) + '\n')
        
        # Write the DataFrame content without including the index or the DataFrame's own headers
        vcf_df = vcf_df[all_headers]  # Reorder DataFrame columns to match the header order
        vcf_df.to_csv(f, sep='\t', index=False, header=False)






def main(vcf_file_path, intervals_file_path, output_file_path):
    # Step 1: Read the VCF file
    vcf_df = read_vcf(vcf_file_path) 

    # Step 2: Read intervals and process them in parallel
    intervals = read_intervals(intervals_file_path)
    processed_alleles = process_vcf_intervals(vcf_df, intervals)

    
    # Step 3: Save the processed allele information
    # Extract the first 9 columns as they contain format-specific information
    vcf_base_df = vcf_df.iloc[:, :9]
    
    # Append processed alleles as new columns
    for individual, alleles in processed_alleles.items():
        # Each individual's alleles are added as a new column
        # Ensure alleles list length matches the number of rows in vcf_base_df
        alleles += ['N/A'] * (len(vcf_base_df) - len(alleles))
        vcf_base_df[individual] = alleles

       
    # Example usage assumes `vcf_base_df` is your base DataFrame and `processed_alleles` is the dictionary with allele data
    print_with_timestamp("Printing the final VCF file")
    save_vcf_with_headers(vcf_base_df, output_file_path)
    print_with_timestamp("Done")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python script.py <vcf_file_path> <intervals_file_path> <output_file_path>")
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
