import re
import sys

# Regular expressions to match the relevant log lines
insert_stats_regex = r"===== Insert Stats =====\n(\d+) messages sent\n(\d+) message bytes sent\n(\d+) Lamport clock \(latencies\)\n(\d+) local AES operations\n(\d+) milliseconds wall clock time\n\{(\d+);(\d+);(\d+)\} nanoseconds \{real;user;system\}\nMem: (\d+) KiB"
extract_stats_regex = r"===== Extract Min Stats =====\n(\d+) messages sent\n(\d+) message bytes sent\n(\d+) Lamport clock \(latencies\)\n(\d+) local AES operations\n(\d+) milliseconds wall clock time\n\{(\d+);(\d+);(\d+)\} nanoseconds \{real;user;system\}\nMem: (\d+) KiB"

# Function to parse insert and extract stats
def parse_logs(log_file):
    with open(log_file, "r") as file:
        log_data = file.read()

    insert_stats = re.findall(insert_stats_regex, log_data)
    extract_stats = re.findall(extract_stats_regex, log_data)

    return insert_stats, extract_stats

# Function to print stats table
def print_stats_table(stats):
    print("Messages Sent | Message Bytes Sent | Lamport Clock | Local AES Operations | Wall Clock Time (ms) | Real Time | User Time | System Time | Memory (KiB)")
    print("-" * 117)
    for stat in stats:
        print("{:14} | {:19} | {:13} | {:21} | {:21} | {:10} | {:9} | {:12} | {:12}".format(*stat))
    print()

if len(sys.argv) != 2:
    print("Usage: python log_parser.py <log_file>")
    sys.exit(1)

log_file = sys.argv[1]

# Parse logs and print stats table
insert_stats, extract_stats = parse_logs(log_file)

print("Insert Stats:")
print_stats_table(insert_stats)

print("Extract Min Stats:")
print_stats_table(extract_stats)

