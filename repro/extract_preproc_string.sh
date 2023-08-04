#!/bin/bash

file_path="$1"

extract_values() {
    while IFS= read -r line; do
        if [[ $line == "Precomputed values used:"* ]]; then
            echo "${line#Precomputed values used: }"
            break
        fi
    done < "$file_path"
}

if [ -z "$file_path" ]; then
    echo "Usage: ./extract_values.sh <file_path>"
else
    extract_values
fi

