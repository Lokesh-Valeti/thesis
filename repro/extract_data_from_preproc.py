import argparse

def parse_file(file_path):
    data = {'heapsize': None, 'is_optimized': None, 'P0': {}, 'P1': {}, 'P2': {}}
    current_party = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("heapsize:"):
                data['heapsize'] = int(line.split(":")[1].strip())
            elif line.startswith("is_optimized:"):
                data['is_optimized'] = int(line.split(":")[1].strip())
            elif line.startswith("===== P"):
                current_party = line.split()[1]
            elif current_party and "milliseconds wall clock time" in line:
                time_parts = line.split()
                data[current_party]['wall_clock_time'] = int(time_parts[0])
            elif current_party and "messages sent" in line:
                data[current_party]['messages_sent'] = int(line.split()[0])
            elif current_party and "message bytes sent" in line:
                data[current_party]['message_bytes_sent'] = int(line.split()[0])
    return data

def main():
    parser = argparse.ArgumentParser(description="Parse file and extract heapsize, is_optimized, wall clock time, messages sent, and message bytes sent for all three parties.")
    parser.add_argument("file_path", type=str, help="Path to the input file")
    args = parser.parse_args()
    
    file_path = args.file_path
    parsed_data = parse_file(file_path)
    
    print(f"{parsed_data['heapsize']}")
    print(f"{parsed_data['is_optimized']}")
    #for party, party_data in parsed_data.items():
     #   if party != 'heapsize' and party != 'is_optimized':
      #      print(f"Party: {party}")
      #      if 'wall_clock_time' in party_data:
       #         print(f"Wall Clock Time: {party_data['wall_clock_time']} ms")
       #     if 'messages_sent' in party_data:
        #        print(f"Messages Sent: {party_data['messages_sent']}")
        #    if 'message_bytes_sent' in party_data:
         #       print(f"Message Bytes Sent: {party_data['message_bytes_sent']}")
         #   print()
   # max_wall_clock_time = 0
   # max_message_bytes_sent = 0

    max_wall_clock_time = 0
    max_message_bytes_sent = 0

    for party, party_data in parsed_data.items():
        if party != 'heapsize' and party != 'is_optimized':
            #print(f"Party: {party}")
            if 'wall_clock_time' in party_data:
                wall_clock_time = party_data['wall_clock_time']
               # print(f"Wall Clock Time: {wall_clock_time} ms")
                max_wall_clock_time = max(max_wall_clock_time, wall_clock_time)
            if 'messages_sent' in party_data:
                messages_sent = party_data['messages_sent']
                #print(f"Messages Sent: {messages_sent}")
            if 'message_bytes_sent' in party_data:
                message_bytes_sent = party_data['message_bytes_sent']
                #print(f"Message Bytes Sent: {message_bytes_sent}")
                max_message_bytes_sent = max(max_message_bytes_sent, message_bytes_sent)
            #print()

    print(f"{max_wall_clock_time}")
    print(f"{max_message_bytes_sent}")


if __name__ == "__main__":
    main()

