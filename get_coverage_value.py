import pandas as pd
import sys

def parse_output_file(output_filename):
    # Initialize an empty list to store the data
    data = []

    # Function to parse the percent identities from the text
    def parse_percent_identities(lines):
        identities = {}
        for i, line in enumerate(lines):
            if 'Percent Identity/Similarity from Consensus Sequence to Reference Genome' in line:
                identities['Percent Identity (Including Gaps)'] = float(lines[i + 1].strip())
                identities['Percent Identity (Excluding Gaps)'] = float(lines[i + 2].strip())
        return identities

    # Open the text file
    with open(output_filename, 'r') as file:
        lines = file.readlines()
        
        # Parse the percent identities
        identities = parse_percent_identities(lines)
        
        # Find the index of the coverage data header
        coverage_start_index = None
        for i, line in enumerate(lines):
            if line.startswith('#rname'):
                coverage_start_index = i + 1
                break

        # Process each subsequent line for coverage data
        if coverage_start_index is not None:
            for line in lines[coverage_start_index:]:
                fields = line.strip().split()
                if len(fields) == 9:
                    try:
                        # Convert fields to the appropriate data types and store in a dictionary
                        record = {
                            'rname': fields[0],
                            'startpos': int(fields[1]),
                            'endpos': int(fields[2]),
                            'numreads': int(fields[3]),
                            'covbases': int(fields[4]),
                            'coverage': float(fields[5]),
                            'meandepth': float(fields[6]),
                            'meanbaseq': float(fields[7]),
                            'meanmapq': float(fields[8])
                        }
                        data.append(record)
                    except ValueError as e:
                        print(f"Error parsing line: {line.strip()}")

    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(data)

    # Add the percent identity values to the DataFrame
    df['Percent Identity (Including Gaps)'] = identities.get('Percent Identity (Including Gaps)', None)
    df['Percent Identity (Excluding Gaps)'] = identities.get('Percent Identity (Excluding Gaps)', None)

    # If there is data in the DataFrame, extract and print the required values
    if not df.empty:
        # Extract the coverage value
        coverage_value = df.loc[0, 'coverage']  # Access the coverage value from the first row

        # Extract the Percent Identity with Gaps
        percent_id_wt_gaps = df.loc[0, 'Percent Identity (Including Gaps)']  # Access the coverage value from the first row

        # Extract the Percent Identity without Gaps
        percent_id_wo_gaps = df.loc[0, 'Percent Identity (Excluding Gaps)']  # Access the coverage value from the first row

        #print(coverage_value)
        #print(percent_id_wt_gaps)
        #print(percent_id_wo_gaps)

        # Create a DataFrame for these values
        summary_df = pd.DataFrame({
            'Coverage': [coverage_value],
            'Percent Identity (Including Gaps)': [percent_id_wt_gaps],
            'Percent Identity (Excluding Gaps)': [percent_id_wo_gaps]
        })

        # Save the summary DataFrame to a CSV file
        summary_csv_filename = f'summary_output_{output_filename}.csv'
        summary_df.to_csv(summary_csv_filename, index=False)

        #print(summary_df)

    return df

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 parser.py <output_filename>")
        sys.exit(1)
    
    output_filename = sys.argv[1]
    parse_output_file(output_filename)
