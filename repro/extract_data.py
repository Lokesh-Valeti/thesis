import re
import sys

def parse_stats(input_string):
    # Regular expressions to extract relevant stats
    heapsize_pattern = r'heapsize: (\d+)'
    insert_stats_pattern = r'===== Insert Stats =====\n([\s\S]*?)\n\n'
    extract_stats_pattern = r'===== Extract Min Stats =====\n([\s\S]*?)(?:\n\n|\Z)'
    optimized_pattern = r'is_optimized: (\d)'

    heapsize = re.search(heapsize_pattern, input_string).group(1)
    insert_stats = re.search(insert_stats_pattern, input_string).group(1)

    # Handling the case when "Extract Min Stats" is not present
    extract_stats_match = re.search(extract_stats_pattern, input_string)
    if extract_stats_match:
        extract_stats = extract_stats_match.group(1)
    else:
        extract_stats = ""

    optimized = re.search(optimized_pattern, input_string).group(1)

    # Extracting insert and extract statistics
    insert_stat_values = re.findall(r'(\d+) messages sent\n(\d+) message bytes sent\n(\d+) Lamport clock \(latencies\)\n(\d+) local AES operations\n(\d+) milliseconds wall clock time', insert_stats)
    extract_stat_values = re.findall(r'(\d+) messages sent\n(\d+) message bytes sent\n(\d+) Lamport clock \(latencies\)\n(\d+) local AES operations\n(\d+) milliseconds wall clock time', extract_stats)

    return {
        "heapsize": heapsize,
        "is_optimized": optimized,
        "insert_stats": insert_stat_values,
        "extract_stats": extract_stat_values
    }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_data.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]

    try:
        with open(filename, "r") as file:
            input_string = file.read()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    parsed_data = parse_stats(input_string)

    insert_stats_wallclock  = 0
    insert_stats_message_bytes = 0
    if parsed_data["insert_stats"]:
        stat = parsed_data["insert_stats"][0]
        insert_stats_wallclock = stat[4]
        insert_stats_message_bytes = stat[1]
   
    extract_stats_wallclock  = 0
    extract_stats_message_bytes = 0

    if parsed_data["extract_stats"]:
        stat = parsed_data["extract_stats"][0]
        extract_stats_wallclock  = stat[4]
        extract_stats_message_bytes = stat[1]

    print(parsed_data["heapsize"])
    print(parsed_data["is_optimized"])
    print(insert_stats_wallclock)
    print(insert_stats_message_bytes)
    print(extract_stats_wallclock)
    print(extract_stats_message_bytes)
