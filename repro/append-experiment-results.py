import sys

# Define a function to extract variables from the data file
def extract_variables_from_file(file_path):
    variables_list = []

    with open(file_path, 'r') as file:
        header = file.readline().strip()

        for line in file:
            columns = line.strip().split()
            ds, op, is_optimized, heapsize, y, pm, y_err = columns

            variables_dict = {
                'DS': ds,
                'OP': op,
                'is_optimized': is_optimized,
                'heapsize': heapsize,
                'y': float(y),  # Convert 'y' to a float to handle averages
                'pm': pm,
                'y_err': float(y_err)  # Convert 'y_err' to a float to handle averages
            }

            variables_list.append(variables_dict)

    return variables_list

# Define a function to update the data based on input and add new rows if needed
def update_data(file_path, ds, op, is_optimized, heapsize, y0):
    extracted_variables = extract_variables_from_file(file_path)

    # Find the row with matching 'DS', 'OP', 'is_optimized', and 'heapsize'
    matching_row = None
    for variables in extracted_variables:
        if (variables['DS'] == ds and
            variables['OP'] == op and
            variables['is_optimized'] == is_optimized and
            variables['heapsize'] == heapsize):
            matching_row = variables
            break

    if matching_row:
        # Update 'y' and 'y_err' based on 'y0'
        matching_row['y'] = (matching_row['y'] + y0) / 2
        matching_row['y_err'] = abs(matching_row['y'] - y0)
    else:
        # If no matching row found, add a new row
        new_row = {
            'DS': ds,
            'OP': op,
            'is_optimized': is_optimized,
            'heapsize': heapsize,
            'y': y0,
            'pm': '±',
            'y_err': 0  # Since it's a new row, we set y_err to 0
        }
        extracted_variables.append(new_row)

    # Write the updated data back to the file
    with open(file_path, 'w') as file:
        header_line = "DS OP is_optimized heapsize y pm y-err\n"
        file.write(header_line)

        for variables in extracted_variables:
            row = f"{variables['DS']} {variables['OP']} {variables['is_optimized']} {variables['heapsize']} {variables['y']} {variables['pm']} {variables['y_err']}\n"
            file.write(row)

# Check if the user provided the required input
#if len(sys.argv) != 6:
 #   print("Usage: python3 script.py <file_path> <DS> <OP> <is_optimized> <heapsize> <y0>")
 #   sys.exit(1)

# Get the input values from the command-line arguments
file_path = sys.argv[1]
ds = sys.argv[2]
op = sys.argv[3]
is_optimized = sys.argv[4]
heapsize = sys.argv[5]
y0 = float(sys.argv[6])  # Convert y0 to a float

# Call the function to update the data based on input
update_data(file_path, ds, op, is_optimized, heapsize, y0)
