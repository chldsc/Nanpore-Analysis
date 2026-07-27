import pysam
import numpy as np
import matplotlib.pyplot as plt

def count_formatter(count):
    if count >= 1_000_000:
        return f"{count/1_000_000:.2f}M"
    elif count >= 1_000:
        return f"{count/1_000:.2f}K"
    else:
        return str(int(count))
    
#calculates match percentage for one read
def match_percentage(counter, read, start_region, end_region, debug = False):
    aligned_bases = 0
    total_length = (end_region - start_region)
    aligned_pairs = read.get_aligned_pairs(matches_only=True, with_seq=True)
    if debug == True:
        print(f"\nRead {counter}")
    for query_pos, ref_pos, base in aligned_pairs:
        if ref_pos is not None and start_region <= ref_pos < end_region:
            #count mutations
            if debug == True:
                print(f"Reference Position: {ref_pos + 1}", f"Base: {base}")
            #lowercase/mismatch is skipped
            if base != 's' and base == base.lower():
                continue
            else:
                aligned_bases += 1
    region_size = end_region - start_region
    final_percent = (aligned_bases/region_size)*100
    if debug == True:
        print(f"Percentage Match: {final_percent}")
    return (aligned_bases/region_size)*100
    
def frequency_pipeline(bamfile, contig, start_region, end_region, total_reads, start_contig, end_contig, debug = False):
    #region_size = end_region - start_region
    start_region = start_region - 1
    bam = pysam.AlignmentFile(bamfile, "rb")

    total_reads_contig = 0
    unmapped_reads = 0
    zero_len_reads = 0
    counter = 1
    final_percentages = []

    # Regions are ANY read that has a read in that region 
    for read in bam.fetch(contig, start_contig, end_contig): #hard coded rn
        total_reads_contig += 1   
        
        if read.is_unmapped:
            unmapped_reads += 1
            continue
        if read.query_length == 0:
            zero_len_reads += 1
            final_percentages.append(0)
            continue
            
        # Calculate match percentage for given read
        read_percent = match_percentage(counter, read, start_region, end_region, debug = debug)
        final_percentages.append(read_percent)
        counter +=1
    not_in_contig = total_reads-total_reads_contig
    for i in range(not_in_contig):
        final_percentages.append(0)
        
    
    print(f"Unmapped Reads: {unmapped_reads} \
            \nZero Length Reads: {zero_len_reads} \
            \nTotal Reads Analyzed: {total_reads} \
            \nReads Absent from Contig: {not_in_contig} \
            \nTotal Reads Contig: {total_reads_contig}")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    _, _, bars = axes[0].hist(final_percentages, color="#1f77b4", edgecolor="black", bins=np.arange(0, 110, 10))
    axes[0].bar_label(bars, color="red", labels=[count_formatter(bar.get_height()) for bar in bars])
    axes[0].set_title(f"Read Match % (Counts) from {total_reads} reads")
    axes[0].set_xlabel("Read Match %")
    axes[0].set_ylabel("Number of Reads")
    axes[0].set_xticks(np.arange(0, 110, 10))

    _, _, bars = axes[1].hist(final_percentages, color="#1f77b4", edgecolor="black",
                              bins=np.arange(0, 110, 10),
                              weights=np.ones(len(final_percentages)) / len(final_percentages))
    axes[1].bar_label(bars, color="red", fmt='{:.2%}')
    axes[1].set_title(f"Read Match % (Proportions) from {total_reads} reads")
    axes[1].set_xlabel("Read Match %")
    axes[1].set_ylabel("Proportion of Reads")
    axes[1].set_xticks(np.arange(0, 110, 10))
    
    plt.tight_layout()

def threshold_processor(threshold, inpath, outpath):
    #generates filtered file and returns number of final reads and prints number filtered out
    inbam = pysam.AlignmentFile(inpath, "rb")
    outbam = pysam.AlignmentFile(outpath, "wb", template=inbam)
    initial_reads = 0
    final_reads = 0
    for read in inbam.fetch(until_eof=True):
        initial_reads += 1
        if read.query_length >= threshold:
            outbam.write(read)
            final_reads += 1
    
    inbam.close()
    outbam.close()
    print(f"Number of Reads Below Threshold: {initial_reads - final_reads} \n\
    Total Reads: {final_reads}")
    return final_reads

def itr_classifier(bamfile, contig, percent_threshold, five_start, five_end, three_start, three_end, total_reads, start_contig, end_contig, debug = False):
    five_start = five_start - 1
    three_start = three_start - 1
    bam = pysam.AlignmentFile(bamfile, "rb")

    total_reads_contig = 0
    unmapped_reads = 0
    zero_len_reads = 0
    counter = 1
    five_prime_percents = []
    three_prime_percents = []
    final_counts = [0,0,0,0]

    # Regions are ANY read that has a read in that region 
    for read in bam.fetch(contig, start_contig, end_contig): #hard coded rn
        total_reads_contig += 1   
        
        if read.is_unmapped:
            unmapped_reads += 1
            continue
        if read.query_length == 0:
            zero_len_reads += 1
            final_counts[0]+=1
            continue
            
        # Calculate match percentage for given read
        five_percent = match_percentage(counter, read, five_start, five_end, debug = debug)
        three_percent = match_percentage(counter, read, three_start, three_end, debug = debug)
        #0 = none, 1= just five prime, 2= just three prime 3= both
        if five_percent < percent_threshold and three_percent < percent_threshold:
            final_counts[0]+=1
        elif five_percent >= percent_threshold and three_percent < percent_threshold:
            final_counts[1]+=1
        elif five_percent < percent_threshold and three_percent >= percent_threshold:
            final_counts[2]+=1
        else:
            final_counts[3]+=1
        counter +=1

    not_in_contig = total_reads-total_reads_contig
    for i in range(not_in_contig):
        final_counts[0]+=1
    '''
    category_labels = {
    0: "None",
    1: "5' Only",
    2: "3' Only",
    3: "Both" }

    counts = {cat: final_categories.count(cat) for cat in category_labels}
    labels= [f"{v}" for v in counts.values()]
    '''
    categories = ["None", "'5' Only", "3' Only", "Both"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    bars = axes[0].bar(categories, final_counts, color="#1f77b4", edgecolor="black")
    axes[0].bar_label(bars, color="red", labels=[count_formatter(bar.get_height()) for bar in bars])
    axes[0].set_title(f"Read Match % (Counts) from {total_reads} reads")
    axes[0].set_xlabel("Read Match %")
    axes[0].set_ylabel("Number of Reads")

    proportions = [count / total_reads for count in final_counts]
    bars = axes[1].bar(categories, proportions, color="#1f77b4", edgecolor="black")
    axes[1].bar_label(bars, color="red", fmt='{:.2%}')
    axes[1].set_title(f"Read Match % (Proportions) from {total_reads} reads")
    axes[1].set_xlabel("Read Match %")
    axes[1].set_ylabel("Proportion of Reads")

def egfp_cassette_by_qpcr(bamfile, contig, contig_coords, percent_threshold, coordinates, total_reads, debug = False):
    #Separate gfp qpcr into yes or no
    #then separate into egfp yes or no
    #then also calculate cassette yes or no 
    #alternates: just egfp, egfp and cassette, neither for gfp and qpcr
    #coordinates format: [(qpcr),(egfp),(cassette)]
    qpcr_start, qcpr_end = coordinates[0][0]-1, coordinates[0][1]
    egfp_start, egfp_end = coordinates[1][0]-1, coordinates[1][1]
    cassette_start, cassette_end = coordinates[2][0]-1, coordinates[2][1]
    start_contig, end_contig = contig_coords[0], contig_coords[1]
    bam = pysam.AlignmentFile(bamfile, "rb")

    total_reads_contig = 0
    yes_qpcr = 0
    yes_qpcr_yes_egfp = 0
    yes_qpcr_no_egfp = 0
    yes_qpcr_yes_cassette = 0
    yes_qpcr_no_cassette = 0
    no_qpcr = 0
    no_qpcr_yes_egfp = 0
    no_qpcr_no_egfp = 0
    no_qpcr_yes_cassette = 0
    no_qpcr_no_cassette = 0   
    counter = 1

    for read in bam.fetch(contig, start_contig, end_contig):
        total_reads_contig += 1   
        
        if read.is_unmapped:
            unmapped_reads += 1                                                                 
            continue
        if read.query_length == 0:
            no_qpcr+=1
            no_qpcr_no_egfp += 1
            no_qpcr_no_egfp += 1
            continue
            
        # Calculate match percentage for given read
        qpcr_percent = match_percentage(counter, read, qpcr_start, qcpr_end, debug = debug)
        #there IS relatively intact qpcr
        if qpcr_percent >= percent_threshold:
            yes_qpcr+=1 

            if match_percentage(counter, read, egfp_start, egfp_end, debug = debug) >= percent_threshold:
                yes_qpcr_yes_egfp +=1
            else: 
                yes_qpcr_no_egfp+=1

            if match_percentage(counter, read, cassette_start, cassette_end, debug = debug) >= percent_threshold:
                yes_qpcr_yes_cassette+=1
            else:
                yes_qpcr_no_cassette+=1

        #there ISNT relatively intact qpcr
        else:
            no_qpcr+=1

            if match_percentage(counter, read, egfp_start, egfp_end, debug = debug) >= percent_threshold:
                no_qpcr_yes_egfp +=1
            else: 
                no_qpcr_no_egfp+=1

            if match_percentage(counter, read, cassette_start, cassette_end, debug = debug) >= percent_threshold:
                no_qpcr_yes_cassette+=1
            else:
                no_qpcr_no_cassette+=1

        counter +=1

    not_in_contig = total_reads-total_reads_contig
    for i in range(not_in_contig):
        no_qpcr += 1
        no_qpcr_no_egfp += 1
        no_qpcr_no_cassette+= 1
    print(f"Total Reads Analyzed: {total_reads}\n\
    Reads with qpcr and EGFP: Count = {yes_qpcr_yes_egfp}\n\
    Reads with qpcr and NO EGFP: Count = {yes_qpcr_no_egfp}\n\
    Reads with qpcr and Cassette: Count = {yes_qpcr_yes_cassette}\n\
    Reads with qpcr and NO Cassette: Count = {yes_qpcr_no_cassette}\n\
        \n\
    Reads with no qpcr and EGFP: Count = {no_qpcr_yes_egfp}\n\
    Reads with no qpcr and NO EGFP: Count = {no_qpcr_no_egfp}\n\
    Reads with no qpcr and Cassette: Count = {no_qpcr_yes_cassette}\n\
    Reads with no qpcr and NO Cassette: Count = {no_qpcr_no_cassette}\n\
    Times reads were looked at:{counter}")
