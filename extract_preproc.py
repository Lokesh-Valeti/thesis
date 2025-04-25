import re
import sys

def extract_preprocessing_strings(filename):
    preprocessing_strings = []

    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'Precomputed values used: (.+)', line)
            if match:
                preprocessing_strings.append(match.group(1))

    for s in preprocessing_strings:
        print(s)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_preprocessing.py <filename>")
    else:
        extract_preprocessing_strings(sys.argv[1])

