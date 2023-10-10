import argparse

overlaps_to_info = {} # connect prefix/suffix to a a list of tuples : [(forward/reverse as a boolean 1:fwd/0:rc, the node id, the unitig's leftmost kmer, prefix/suffix as a boolean 1:pref/0:suff)] we use a list beauce k-1 mers can happen several times
kmers_to_edge = {} # connect a kmer (the first kmer read in the unitig) to the unitigs edges info



def reverse_complement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in reversed(dna)])

def is_canonical(kmer):
    rc = reverse_complement(kmer)
    if kmer < rc:
        return (True, kmer)
    else:
        return (False, rc)



def parse_fasta_buffered(file_path, k, left_overlap_list, right_overlap_list, buffer_size=1000):
    """Parse FASTA file and fills to global dictionaries + returns sorted lists of k-1 mers to find connections between unitigs

    Args:
    - file_path (str): Path to the FASTA file.
    - k (int): kmer size
    - left_overlap_list (list): an overlap list
    - right_overlap_list (list): an overlap list
    - buffer_size (int): Number of lines to read at once.

    Returns:
    two sorted lists of k-1 mers
    """
    header = None
    buffer = []

    with open(file_path, 'r') as f:
        node_id = 0
        while True:
            # Read lines into buffer
            while len(buffer) < buffer_size:
                line = f.readline()
                if not line:
                    break
                buffer.append(line.rstrip())

            if not buffer:
                break

            while buffer:
                line = buffer.pop(0)
                if line.startswith(">"):  # Header line
                    header = line
                else:
                    left_kmer = line[:k]
                    right_kmer = line[len(line)-k:]
                    is_cl, left_canonical = is_canonical(left_kmer[:-1])
                    is_cr, right_canonical = is_canonical(right_kmer[1:])

                    several = overlaps_to_info.get(left_canonical)
                    if several is not None and len(several) > 0 :
                        overlaps_to_info[left_canonical].append((True, is_cl, node_id, left_kmer)) # True: prefix/False: suffix, True: canon/False: rc, node_id (int)
                    else:
                        overlaps_to_info.setdefault(left_canonical, [(True, is_cl, node_id, left_kmer)])
                    # insert_overlaps(left_overlap_list, right_overlap_list, left_canonical, is_cl, True) #put k-1 mer in list

                    several = overlaps_to_info.get(right_canonical)
                    if several is not None and len(several) > 0:
                        overlaps_to_info[right_canonical].append((False, is_cr, node_id, left_kmer))
                    else:
                        overlaps_to_info.setdefault(right_canonical,[(False, is_cr, node_id, left_kmer)])
                    # insert_overlaps(left_overlap_list, right_overlap_list, right_canonical, is_cr, False) #put k-1 mer in list
                    # print(line)
                    # print(overlaps_to_info)
                    kmers_to_edge[left_kmer] = ""
                    node_id += 1
                    


            # Handle the remaining sequence in the buffer for the current header
                # if header and buffer:
                #     line = buffer
                #     buffer = []
                #     left_kmer = line[:k]
                #     right_kmer = line[len(line)-k:]
                #     is_cl, left_canonical = is_canonical(left_kmer[:-1])
                #     is_cr, right_canonical = is_canonical(right_kmer[1:])
                #     if left_kmer == "AAAAACTATGGATGTCAAAAAACCCGGCAA" or right_kmer == "AAAAACTATGGATGTCAAAAAACCCGGCAA":
                #         print("ID", node_id)
                #     overlaps_to_info[left_canonical] = (0, node_id, left_kmer)
                #     insert_overlaps(left_overlap_list, right_overlap_list, left_canonical, is_cl) #put k-1 mer in list
                #     overlaps_to_info[right_canonical] = (1, node_id, left_kmer)
                #     insert_overlaps(left_overlap_list, right_overlap_list, right_canonical, is_cr) #put k-1 mer in list
                #     kmers_to_edge[left_kmer] = ""
                #     node_id += 1
    f.close()

    return left_overlap_list, right_overlap_list


    
def find_connections2():
    for overlap in overlaps_to_info.keys():
            node1 = overlaps_to_info[overlap][0]
            for node2 in overlaps_to_info[overlap][1:]:
                if node1[0] is True: #prefix
                    if node1[1] is True: #canonical prefix
                        if (node2[0] is True) and (node2[1] is False): #reverse prefix: possible link
                            kmers_to_edge[node1[3]] += "L:+:"+ str(node2[2]) + ":- "
                            if (node2[2] != node1[2]): #avoid self connected edge case
                                kmers_to_edge[node2[3]] += "L:-:"+ str(node1[2]) + ":+ "
                        elif (node2[0] is False) and (node2[1] is True): # canon suffix: possible link
                            kmers_to_edge[node1[3]] += "L:+:"+ str(node2[2]) + ":+ "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:+:"+ str(node1[2]) + ":+ "
                    else: #reverse prefix
                        if (node2[0] is True) and (node2[1] is True): #canon prefix: possible link
                            kmers_to_edge[node1[3]] += "L:-:"+ str(node2[2]) + ":+ "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:+:"+ str(node1[2]) + ":- "
                        elif (node2[0] is False) and (node2[1] is False): # reverse suffix: possible link
                            kmers_to_edge[node1[3]] += "L:-:"+ str(node2[2]) + ":- "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:-:"+ str(node1[2]) + ":- "
                else: #suffix
                    if node1[1] is False: #canonical suffix
                        if (node2[0] is True) and (node2[1] is True): #canon prefix: possible link
                            kmers_to_edge[node1[3]] += "L:+:"+ str(node2[2]) + ":+ "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:+:"+ str(node1[2]) + ":+ "
                        elif (node2[0] is False) and (node2[1] is False): # reverse suffix: possible link
                            kmers_to_edge[node1[3]] += "L:+:"+ str(node2[2]) + ":- "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:-:"+ str(node1[2]) + ":+ "

                    else: #reverse suffix
                        if (node2[0] is True) and (node2[1] is False): #reverse prefix: possible link
                            kmers_to_edge[node1[3]] += "L:-:"+ str(node2[2]) + ":- "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:-:"+ str(node1[2]) + ":- "
                        elif (node2[0] is False) and (node2[1] is True): # canon suffix: possible link
                            kmers_to_edge[node1[3]] += "L:-:"+ str(node2[2]) + ":+ "
                            if (node2[2] != node1[2]):
                                kmers_to_edge[node2[3]] += "L:+:"+ str(node1[2]) + ":- "


               





def write_bcalm(file_path, output_path, k, buffer_size=1000):
    """Does another pass on the input file and rewrites the headers to get the connections
    Args:
    - file_path (str): input file
    - output_path (str): output file
    - k (int): k size
    """
    node_id = 0
    header = None
    buffer = []
    with open(output_path, 'w') as o:
        with open(file_path, 'r') as f:
            while True:
                # Read lines into buffer
                while len(buffer) < buffer_size:
                    line = f.readline()
                    if not line:
                        break
                    buffer.append(line.rstrip())

                if not buffer:
                    break

                while buffer:
                    header = buffer.pop(0)
                    elements_header = header.split(">") #gets the color from ggcat
                    sequence = buffer.pop(0)
                    pref_kmer = sequence[:k]
                    o.write(">" + str(node_id) + " LN:i:0 KC:i:0 km:f:0.0\t" + kmers_to_edge[pref_kmer] + "\t" + elements_header[1] +"\n"+ sequence + "\n")
                    node_id += 1
                
                if header and buffer:
                    o.write(buffer)
                    node_id += 1
    
    o.close()
    f.close()



def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Create BCALM from GGCAT output")

    # Add arguments
    parser.add_argument('-k', '--kmersize', type=int, help='k-mer value')
    parser.add_argument('-f', '--input', type=str, help='input ggcat file')
    parser.add_argument('-o', '--output', type=str, help='output file for bcalm')

    # Parse arguments
    args = parser.parse_args()

    k = args.kmersize
    input_file = args.input
    output_file = args.output

    left_overlap_list = list()
    right_overlap_list = list()
    headers = list()
    left_overlap_list, right_overlap_list = parse_fasta_buffered(input_file, k, left_overlap_list, right_overlap_list)
    find_connections2()
    # find_connections(left_overlap_list, right_overlap_list)
    write_bcalm(input_file, output_file, k)

if __name__ == "__main__":
    main()
    
