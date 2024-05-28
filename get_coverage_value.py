# Initialize an empty list to store the data
data = []

# Open the text file
with open('output.text', 'r') as file:
    # Skip the header line
    header = file.readline().strip().split()

    # Process each subsequent line
    for line in file:
        # Split the line by tabs or spaces
        fields = line.strip().split()
        
        # Verify that the line has the expected number of fields
        if len(fields) == 9:
            try:
                # Convert fields to the appropriate data types and store in a dictionary
                record = {
                    'rname': fields[0],
                    'startpos': int(fields[1]),
                    'endpos': int(fields[2]),
                    'numreads': int(fields[3]),
                    'covbases': int(fields[4]),
                    'coverage': float(fields[5]),  # Assuming coverage might be a float
                    'meandepth': float(fields[6]),
                    'meanbaseq': float(fields[7]),
                    'meanmapq': float(fields[8])
                }
                data.append(record)
            except ValueError as e:
                print("")
        else:
            print("")

# Create a DataFrame from the list of dictionaries
df = pd.DataFrame(data)

# Extract the coverage value
coverage_value = df.loc[0, 'coverage']  # Access the coverage value from the first row

# Save the DataFrame to a CSV file
df.to_csv('output_table.csv', index=False)

# Print the DataFrame
print(df)
print(coverage_value)