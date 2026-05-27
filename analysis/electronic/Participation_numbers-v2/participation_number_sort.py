# Import necessary module

# Define a function to reorder lines by the third column
def reorder_lines_by_third_column(input_file, output_file):
    try:
        # Read lines from the input file
        with open(input_file, 'r') as file:
            lines = file.readlines()

        # Parse the lines into a list of tuples with the third column as a key
        lines_parsed = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:  # Ensure there are at least three columns
                lines_parsed.append((line, float(parts[2])))

        # Sort lines by the third column in descending order
        lines_sorted = sorted(lines_parsed, key=lambda x: x[1], reverse=True)

        # Write the sorted lines back to the output file
        with open(output_file, 'w') as file:
            for line, _ in lines_sorted:
                file.write(line)

        print(f"Reordered lines written to {output_file}.")

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
# Input: "Participation_number.dat", Output: "output.txt"
input_file = "Participation_number.dat"
output_file = "Participation_number_sorted.dat"
reorder_lines_by_third_column(input_file, output_file)

