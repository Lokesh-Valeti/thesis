import sys
import re

def extract_binary_search_data(output):
    bytes_sent_pattern = r"BINARY SEARCH =====\n\d+ messages sent\n(\d+) message bytes sent"
    wall_clock_pattern = r"BINARY SEARCH =====\n\d+ messages sent\n\d+ message bytes sent\n\d+ Lamport clock \(latencies\)\n\d+ local AES operations\n(\d+) milliseconds wall clock time"

    bytes_sent_matches = re.findall(bytes_sent_pattern, output)
    wall_clock_matches = re.findall(wall_clock_pattern, output)

    bytes_sent = [int(match) for match in bytes_sent_matches]
    wall_clock_time = [int(match) for match in wall_clock_matches]

    return bytes_sent, wall_clock_time

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    try:
        with open(input_file, 'r') as f:
            output = f.read()
    except FileNotFoundError:
        print(f"File not found: {input_file}")
        sys.exit(1)

    bytes_sent, wall_clock_time = extract_binary_search_data(output)

    max_bytes_sent = max(bytes_sent)
    max_wall_clock_time = max(wall_clock_time)

    print(f"{max_bytes_sent}")
    print(f"{max_wall_clock_time}")

