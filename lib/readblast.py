# for multi_padlock_design
# Xiaoyan, 2017

import os
from lib import createoutput
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

notmapped = []
funmap = []

dna_IMM1 = {
    "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0), "GG/CT": (3.3, 10.4), "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3), "GT/TG": (4.1, 9.5), "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2), "TT/AG": (-1.3, -5.3), "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3), "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2), "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2), "AT/TC": (-1.2, -6.2), "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1), "GC/CT": (2.3, 5.4), "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7), "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6), "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9), "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9), "TC/AC": (6.1, 16.4), "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8), "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
    "AI/TC": (-8.9, -25.5), "TI/AC": (-5.9, -17.4), "AC/TI": (-8.8, -25.4),
    "TC/AI": (-4.9, -13.9), "CI/GC": (-5.4, -13.7), "GI/CC": (-6.8, -19.1),
    "CC/GI": (-8.3, -23.8), "GC/CI": (-5.0, -12.6),
    "AI/TA": (-8.3, -25.0), "TI/AA": (-3.4, -11.2), "AA/TI": (-0.7, -2.6),
    "TA/AI": (-1.3, -4.6), "CI/GA": (2.6, 8.9), "GI/CA": (-7.8, -21.1),
    "CA/GI": (-7.0, -20.0), "GA/CI": (-7.6, -20.2),
    "AI/TT": (0.49, -0.7), "TI/AT": (-6.5, -22.0), "AT/TI": (-5.6, -18.7),
    "TT/AI": (-0.8, -4.3), "CI/GT": (-1.0, -2.4), "GI/CT": (-3.5, -10.6),
    "CT/GI": (0.1, -1.0), "GT/CI": (-4.3, -12.1),
    "AI/TG": (-4.9, -15.8), "TI/AG": (-1.9, -8.5), "AG/TI": (0.1, -1.8),
    "TG/AI": (1.0, 1.0), "CI/GG": (7.1, 21.3), "GI/CG": (-1.1, -3.2),
    "CG/GI": (5.8, 16.9), "GG/CI": (-7.6, -22.0),
    "AI/TI": (-3.3, -11.9), "TI/AI": (0.1, -2.3), "CI/GI": (1.3, 3.0),
    "GI/CI": (-0.5, -1.3),
    
    "AA/AA": (3.8, 6.3), "AA/AC": (3.9, 6.3), "AA/AG": (3.8, 6.3),
    "AA/CA": (3.9, 6.6), "AA/CC": (4.1, 6.6), "AA/CG": (4.4, 6.6),
    "AA/GA": (3.3, 4.7), "AA/GC": (3.4, 4.7), "AA/GG": (3.3, 4.7),
    "AC/AA": (4.3, 7.5), "AC/AC": (4.4, 7.5), "AC/AT": (4.9, 7.5),
    "AC/CA": (4.5, 7.8), "AC/CC": (4.7, 7.8), "AC/CT": (4.6, 7.8),
    "AC/GA": (3.8, 5.9), "AC/GC": (3.9, 5.9), "AC/GT": (4.3, 5.9),
    "AG/AA": (4.1, 7.0), "AG/AG": (4.1, 7.0), "AG/AT": (4.7, 7.0),
    "AG/CA": (4.3, 7.3), "AG/CG": (4.8, 7.3), "AG/CT": (4.4, 7.3),
    "AG/GA": (3.6, 5.4), "AG/GG": (3.6, 5.4), "AG/GT": (4.1, 5.4),
    "AT/AC": (4.6, 6.6), "AT/AG": (4.5, 6.6), "AT/AT": (5.1, 6.6),
    "AT/CC": (4.9, 6.9), "AT/CG": (5.2, 6.9), "AT/CT": (4.8, 6.9),
    "AT/GC": (4.2, 5.0), "AT/GG": (4.1, 5.0), "AT/GT": (4.6, 5.0),
    "CA/AA": (4.0, 6.3), "CA/AC": (4.1, 6.3), "CA/AG": (4.0, 6.3),
    "CA/CA": (4.1, 6.6), "CA/CC": (4.3, 6.6), "CA/CG": (4.6, 6.6),
    "CA/TA": (2.0, -1.8), "CA/TC": (1.6, -1.8), "CA/TG": (2.0, -1.8),
    "CC/AA": (4.5, 7.5), "CC/AC": (4.6, 7.5), "CC/AT": (5.1, 7.5),
    "CC/CA": (4.7, 7.8), "CC/CC": (4.9, 7.8), "CC/CT": (4.8, 7.8),
    "CC/TA": (2.6, -0.6), "CC/TC": (2.2, -0.6), "CC/TT": (2.0, -0.6),
    "CG/AA": (4.3, 7.0), "CG/AG": (4.3, 7.0), "CG/AT": (4.9, 7.0),
    "CG/CA": (4.5, 7.3), "CG/CG": (5.0, 7.3), "CG/CT": (4.6, 7.3),
    "CG/TA": (2.4, -1.1), "CG/TG": (2.4, -1.1), "CG/TT": (1.8, -1.1),
    "CT/AC": (4.2, 6.6), "CT/AG": (4.1, 6.6), "CT/AT": (4.7, 6.6),
    "CT/CC": (4.5, 6.9), "CT/CG": (4.8, 6.9), "CT/CT": (4.4, 6.9),
    "CT/TC": (1.8, -1.5), "CT/TG": (2.2, -1.5), "CT/TT": (1.6, -1.5),
    "GA/AA": (3.9, 6.3), "GA/AC": (4.0, 6.3), "GA/AG": (3.9, 6.3),
    "GA/GA": (3.4, 4.7), "GA/GC": (3.5, 4.7), "GA/GG": (3.4, 4.7),
    "GA/TA": (1.9, -1.8), "GA/TC": (1.5, -1.8), "GA/TG": (1.9, -1.8),
    "GC/AA": (4.6, 7.5), "GC/AC": (4.7, 7.5), "GC/AT": (5.2, 7.5),
    "GC/GA": (4.1, 5.9), "GC/GC": (4.2, 5.9), "GC/GT": (4.6, 5.9),
    "GC/TA": (2.7, -0.6), "GC/TC": (2.3, -0.6), "GC/TT": (2.1, -0.6),
    "GG/AA": (4.1, 7.0), "GG/AG": (4.1, 7.0), "GG/AT": (4.7, 7.0),
    "GG/GA": (3.6, 5.4), "GG/GG": (3.6, 5.4), "GG/GT": (4.1, 5.4),
    "GG/TA": (2.2, -1.1), "GG/TG": (2.2, -1.1), "GG/TT": (1.6, -1.1),
    "GT/AC": (4.6, 6.6), "GT/AG": (4.5, 6.6), "GT/AT": (5.1, 6.6),
    "GT/GC": (4.2, 5.0), "GT/GG": (4.1, 5.0), "GT/GT": (4.6, 5.0),
    "GT/TC": (2.2, -1.5), "GT/TG": (2.6, -1.5), "GT/TT": (2.0, -1.5),
    "TA/CA": (4.6, 6.6), "TA/CC": (4.8, 6.6), "TA/CG": (5.1, 6.6),
    "TA/GA": (4.0, 4.7), "TA/GC": (4.1, 4.7), "TA/GG": (4.0, 4.7),
    "TA/TA": (2.5, -1.8), "TA/TC": (2.1, -1.8), "TA/TG": (2.5, -1.8),
    "TC/CA": (4.6, 7.8), "TC/CC": (4.8, 7.8), "TC/CT": (4.7, 7.8),
    "TC/GA": (3.9, 5.9), "TC/GC": (4.0, 5.9), "TC/GT": (4.4, 5.9),
    "TC/TA": (2.5, -0.6), "TC/TC": (2.1, -0.6), "TC/TT": (1.9, -0.6),
    "TG/CA": (4.9, 7.3), "TG/CG": (5.4, 7.3), "TG/CT": (5.0, 7.3),
    "TG/GA": (4.2, 5.4), "TG/GG": (4.2, 5.4), "TG/GT": (4.7, 5.4),
    "TG/TA": (2.8, -1.1), "TG/TG": (2.8, -1.1), "TG/TT": (2.2, -1.1),
    "TT/CC": (4.3, 6.9), "TT/CG": (4.6, 6.9), "TT/CT": (4.2, 6.9),
    "TT/GC": (3.6, 5.0), "TT/GG": (3.5, 5.0), "TT/GT": (4.0, 5.0),
    "TT/TC": (1.6, -1.5), "TT/TG": (2.0, -1.5), "TT/TT": (1.4, -1.5),
    }

def calc_tm_NN(seq, cseq=None, dnac1=60, Na=75, K=75, Tris=20, Mg=10, dNTPs=0, fmd=20, strict=True):
    if cseq is None:
        cseq = seq
    tm = mt.chem_correction(
        mt.Tm_NN(
            seq=Seq(seq),
            c_seq=Seq(cseq).complement(),
            dnac1=dnac1,
            Na=Na,
            K=K,
            Tris=Tris,
            Mg=Mg,
            dNTPs=dNTPs,
            imm_table=dna_IMM1,
            tmm_table=dna_TMM1,
            strict=strict)
        ,fmd=fmd)
    return tm

def has_gap_or_mismatch(query, subject, ligation_site, start_pos, buffer=3):
    # Split the sequence into left and right arms
    query_left, query_right = split_arms(query, ligation_site, start_pos)
    subject_left, subject_right = split_arms(subject, ligation_site, start_pos)
    # Check for gaps or mismatches within buffer of the ligation site
    # For the left arm, check the right end of the sequence
    if len(query_left) < buffer or len(query_right) < buffer:
        return True
    for i in range(-buffer, 1):
        if query_left[i] != subject_left[i]:
            return True

    # For the right arm, check the left end of the sequence
    for i in range(buffer):
        if query_right[i] != subject_right[i]:
            return True

    return False

def split_arms(sequence, ligation_site, start_pos):
    split_point = ligation_site - start_pos
    left_arm = sequence[:split_point]
    right_arm = sequence[split_point:]
    return left_arm, right_arm

def fill_gaps(query, subject):
    """
    Fills gaps in query or subject sequence with corresponding base from the other sequence.

    Args:
        query (str): The query sequence containing possible gaps.
        subject (str): The subject sequence containing possible gaps.

    Returns:
        tuple: Filled query and subject sequences.
    """
    filled_query = []
    filled_subject = []
    
    # Ensure the sequences are the same length
    if len(query) != len(subject):
        raise ValueError("Query and Subject sequences must be of equal length")

    # Iterate through both sequences to fill gaps
    for q_base, s_base in zip(query, subject):
        if q_base == '-' and s_base != '-':
            filled_query.append(s_base)
            filled_subject.append(s_base)
        elif s_base == '-' and q_base != '-':
            filled_query.append(q_base)
            filled_subject.append(q_base)
        else:
            filled_query.append(q_base)
            filled_subject.append(s_base)

    return ''.join(filled_query), ''.join(filled_subject)

def readblastout(file, armlength, variants, specificity_by_tm=False):
    """ Read the results from blast
    
    Args:
        file: str, path to the blast output file
        armlength: int, length of the arm
        variants: list, list of variants to check for
        
    Returns:
        specific: bool, True if the sequence is specific, False if not
    """
    global funmap
    global notmapped
    specific = True
    mappable = False
    with open(file, "r") as fh:
        noblast = True
        for line in fh:
            noblast = False
            if specific:
                if not ("XR" in line or "XM" in line):  # skip all predicted transcripts
                    if "|" in line:
                        columns = line.split("|")
                        if len(columns) <= 3:
                            hit = columns[1]#.split(".", 1)[0]
                            scores = columns[-1].split(",")
                        else:
                            hit = columns[3]#.split(".", 1)[0]
                            scores = columns[-1].split(",")
                    else:
                        columns = line.split(",")
                        hit = columns[1]#.split(".", 1)[0]
                        scores = columns[2:]
                    if len(scores):
                        # remove the first empty column (somehow appears in some db and blast versions)
                        if "" in scores:
                            scores.remove("")
                        # Filter hits based on Tm of each arm and mismatches near the ligation site
                        ligation_site = armlength + 1
                        #query is the padlock sequence
                        query_seq = columns[12]
                        #subject is the transcript database sequence
                        subject_seq = columns[13].split('\n')[0]
                        query_start = int(columns[6])
                        if specificity_by_tm:
                            # First check if the sequences have gaps or mismatches near the ligation site
                            if not has_gap_or_mismatch(query_seq, subject_seq, ligation_site, query_start):
                                # Then split the sequences into left and right arms using actual start
                                # position of matched sequence
                                query_left, query_right = split_arms(query_seq, ligation_site, query_start)
                                subject_left, subject_right = split_arms(subject_seq, ligation_site, query_start)
                                query_left = query_left.strip()
                                subject_left = subject_left.strip()
                                query_right = query_right.strip()
                                subject_right = subject_right.strip()
                                if query_left:
                                    # Next fill gaps in the sequences. This is being generous
                                    # so that we err on the side of caution with potential non-specific hits
                                    query_left, subject_left = fill_gaps(query_left, subject_left)
                                    # Then check Tm of each arm
                                    tm_left = calc_tm_NN(seq=query_left, cseq=subject_left)
                                    tm_no_mismatch_left = calc_tm_NN(seq=query_left)
                                else:
                                    tm_left = None
                                if query_right:
                                    query_right, subject_right = fill_gaps(query_right, subject_right)
                                    tm_right = calc_tm_NN(seq=query_right, cseq=subject_right)
                                    tm_no_mismatch_right = calc_tm_NN(seq=query_right)
                                else:
                                    tm_right = None
                                # If both arms have Tm > 37, it's a valid probe
                                if (tm_left < 37) or (tm_right < 37):
                                    print("Invalid probe")
                                else:
                                    print("Valid probe")
                                    # First check if variants are provided
                                    if len(variants):
                                        # If they are and the hit is not in them, this is a non specific hit
                                        if hit not in variants:
                                            if isinstance(variants, list):
                                                with open(os.path.join(os.path.dirname(file), 'homology.txt'), 'a') as fsimilar:
                                                        fsimilar.write('%s,%s\n' % (hit, variants[0]))
                                            else:
                                                with open(os.path.join(os.path.dirname(file), 'homology.txt'), 'a') as fsimilar:
                                                        fsimilar.write('%s,%s\n' % (hit, variants))
                                            specific = False
                                        # Otherwise, the hit is specific
                                        else:
                                            # And if it's a perfect match mark it as mappable
                                            if float(scores[0]) == 100 and int(scores[1]) == 2*armlength:
                                                mappable = True
                                    # If no variants are provided, check if the hit is the same as the input sequence
                                    else:
                                        import warnings
                                        warnings.warn('No gene variants searched for due to providing fasta input, \
                                                    only checking if blast hits are the same as the input sequence')
                                        if hit not in file.split('/')[-1]:
                                            specific = False
                                        else:
                                            if float(scores[0]) == 100 and int(scores[1]) == 2*armlength:
                                                mappable = True
                        else:
                            # Use original filtering logic if specificity_by_tm is False
                            if (
                                armlength < int(scores[1])
                                and float(scores[0]) > 80
                                and int(scores[4]) < armlength - 4
                                and int(scores[5]) > armlength + 5
                            ):
                                print("Valid probe")
                                # Proceed with specificity check based on variants
                                if len(variants):
                                    if hit not in variants:
                                        if isinstance(variants, list):
                                            with open(os.path.join(os.path.dirname(file), 'homology.txt'), 'a') as fsimilar:
                                                fsimilar.write('%s,%s\n' % (hit, variants[0]))
                                        else:
                                            with open(os.path.join(os.path.dirname(file), 'homology.txt'), 'a') as fsimilar:
                                                fsimilar.write('%s,%s\n' % (hit, variants))
                                        specific = False
                                    else:
                                        if float(scores[0]) == 100 and int(scores[1]) == 2 * armlength:
                                            mappable = True
                                else:
                                    import warnings
                                    warnings.warn('No gene variants searched for due to providing fasta input, \
                                                only checking if blast hits are the same as the input sequence')
                                    if hit not in file.split('/')[-1]:
                                        specific = False
                                    else:
                                        if float(scores[0]) == 100 and int(scores[1]) == 2 * armlength:
                                            mappable = True
                            else:
                                print("Invalid probe")

                        print(f"Tm left:             {tm_left:.2f} C  Tm right:             {tm_right:.2f} C")
                        print(f"Tm no mismatch left: {tm_no_mismatch_left:.2f} C  Tm no mismatch right: {tm_no_mismatch_right:.2f} C")
                        print("Q left:", query_left, "Q right:", query_right)
                        print("S left:", subject_left, "S right:", subject_right)
                        print(f"Transcript ID: {columns[1]}")
                        print(f"E-value: {columns[10]} \n")
            # unmappable sequences will later be removed from final results
            if not mappable:
                with open(file[0:-10] + '.fasta', 'r') as f:
                    seq = f.readlines()
                    funmap.write('Could not map sequence in ' + file[:-10] + '!\n')
                    funmap.write(seq[1])
                    notmapped.append(int(file[:-10].split('_')[-1])-1)

    # if no blast results returned, ignore the sequence (can be due to error in blast+ or due to no similar sequence)
    if noblast:
        specific = False
    return specific


def getcandidates(listSiteChopped, headers, dirnames, armlength, accession, specificity_by_tm):
    """ Get the candidates for the probes

    Args:
        listSiteChopped (list): list of sites
        headers (list): list of headers
        dirnames (list): list of directories
        armlength (int): arm length
        accession (list): list of accession numbers

    Returns:
        siteCandidates (list): list of candidates
        notMapped (list): list of not mapped sites
    """  
    global notmapped
    global funmap
    siteCandidates = []
    notMapped = []
    t = dirnames[1].split("TempFolder")[1]
    with open(os.path.join(dirnames[0], "2.Unmappable_" + t + ".txt"), "w") as funmap:
        for i, sites in enumerate(listSiteChopped):
            funmap.write("%s\n" % headers[i])

            fname = createoutput.blastinfilename(dirnames[1], headers[i])

            import config
            fastadir = (config.fastadir_mouse, config.fastadir_human)
            with open(fastadir[0] + "/" + "mouse" + ".allheaders.txt", "r") as f:
                Headers = [line.rstrip("\n") for line in f]

            #take file and split to find gene name
            gene_name = fname.split("/")[-1].split(".")[0]
            print("Gene name:", gene_name)
            first_match_line = next((line for line in Headers if gene_name in line), None)

            if first_match_line:
                # Step 2: Extract the gene ID
                gene_id_start = first_match_line.find("gene:") + len("gene:")
                gene_id_end = first_match_line.find(" ", gene_id_start)
                gene_id = first_match_line[gene_id_start:gene_id_end]

                # Step 3: Search for all instances of this gene ID
                matching_lines = [line for line in Headers if gene_id in line]

                # Step 4: Extract the required substrings
                variants = [line[1:line.find(" ")] for line in matching_lines]

                print("Gene ID:", gene_id)
                print("Variants:", variants)
            else:
                print("No match found for the gene_name.")
                variants = []

            notmapped = []
            blast_bw = []
            for j, target in enumerate(sites):
                fblast = fname + "_" + str(j + 1) + "_blast.txt"
                blast_bw.append(readblastout(fblast, armlength, variants, specificity_by_tm))

            # find sequences that are specific enough
            idxspecific = np.nonzero(blast_bw)[0]
            tempCandidates = np.array(sites)
            sitespecific = tempCandidates[idxspecific]

            # write unmappable sites
            notmapped = tempCandidates[notmapped]
            recorded = False
            temp = None
            for j, nomap in enumerate(notmapped):
                if j == 0:
                    funmap.write(
                        "\nUnmapped sequence starting position(s):\n%d" % (nomap + 1)
                    )
                    recorded = True
                else:
                    if temp is not None and nomap == temp + 1:
                        funmap.write("-")
                        recorded = False
                        if j == len(notmapped) - 1:
                            funmap.write("%d" % (nomap + 1))
                    else:
                        if recorded:
                            funmap.write(",%d" % (nomap + 1))
                        else:
                            funmap.write("%d,%d" % (temp + 1, nomap + 1))
                            recorded = True
                temp = nomap
            notMapped.append(notmapped)

            # continuous regions
            if len(sitespecific):
                idxPairStart = np.nonzero(sitespecific[1:] - sitespecific[0:-1] != 1)[0]
                if len(idxPairStart) == 0 and len(
                    sitespecific
                ):  # only one continuous region exists
                    idxPairStart = np.array([0])
                    idxPairEnd = np.array([len(sitespecific) - 1])
                else:
                    if idxspecific[0] == 0:
                        idxPairStart = np.append(
                            idxspecific[0], np.add(idxPairStart, 1)
                        )
                    else:
                        idxPairStart = np.add(idxPairStart, 1)

                    idxPairEnd = np.zeros(len(idxPairStart), dtype=int)
                    for j in range(len(idxPairStart) - 1):
                        p = idxPairStart[j]
                        c = 0
                        while sitespecific[p + c + 1] == sitespecific[p + c] + 1:
                            c += 1
                        idxPairEnd[j] = p + c

                sitePairStart = sitespecific[idxPairStart]
                sitePairEnd = sitespecific[idxPairEnd]
                sitePairEnd[-1] = sitespecific[-1]
                siteCandidates.append(np.vstack((sitePairStart, sitePairEnd)))

            else:  # no usable fragment
                siteCandidates.append(np.zeros((2, 0)))
            funmap.write("\n\n")
    return siteCandidates, notMapped
