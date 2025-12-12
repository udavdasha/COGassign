"""
Module for working with different software and databases. Currently used
software blocks:
    (1) Misc methods
    (2) TMHMM
    (3) Pfam
    (4) Other HMMer files
    (5) COG database
    (6) Genbank assembly and Prodigal
    (7) BLAST (from <get_operon.py>)
------- Version: 3.6
"""
import os, sys, re
import math

#------------------------------------------------------------------------------
#                            (1) Misc methods
#------------------------------------------------------------------------------
def read_plain_list(list_filename):
    """
    Method reads plain list inside the <list_filename> and returns it as a dictionary.
    """
    result = dict()
    list_file = open(list_filename)
    for string in list_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        result[string] = True
    list_file.close()
    return result

def read_ordered_list(list_filename):
    """
    Method reads a plain list inside the <list_filename> and returns it as a list,
    i.e. with the preserved order.
    """
    result = list()
    list_file = open(list_filename)
    for string in list_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        result.append(string)
    list_file.close()
    return result

def read_two_column_assignment(input_filename, report_duplicate = True, allow_multiple_columns = False, as_dict = False):
    """
    Method reads two-column file, where first column is a group label (or query, in case of BLAST/HMMer searchs)
    and the second one is assigned to it somehow. If <allow_multiple_columns>, script will ignore the fact that more
    than 2 columns could exist in file.

    Returns tuple of (assignment, second_column), where <assignment> is a dictionary, where each value from the
    first column is assigned to a list of values to which if refferes in the second column (or a single one, if <as_dict> == True).
   
    """
    assignment = dict()
    second_column = dict()
    duplicated_n = 0
    input_file = open(input_filename)
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        if not allow_multiple_columns and (len(fields) != 2):
            print ("FATAL ERROR: more than two fields detected in line: '%s'" % string)
            sys.exit()
        query = fields[0]
        target = fields[1]
        if as_dict == True:
            if not query in assignment:
                assignment[query] = target
            else:
                print("Duplicated assignment: '%s' already assigned to '%s', assignment '%s' skipped" % (query, assignment[query], target))
        else:
            if not query in assignment:
                assignment[query] = list()
            if len(assignment[query]) != 0: # More than one assignment of query to target
                duplicated_n += 1            
            assignment[query].append(target)
           
        second_column[target] = False
    input_file.close()
    if report_duplicate:
        print ("Number of duplicated assignments detected: %i" % duplicated_n)
    return (assignment, second_column)

def vary_motif(motif):
    motif_code = motif.split("_")[0]
    coded_motifs = list()
    coded_motifs.append(motif_code + "_")
    motif_seq = motif.split("_")[1]
    variant_reading = False
    v = 0
    new_motifs = list()
    for l in motif_seq:
        if l == "[":
            variant_reading = True
            v = 0
            new_motifs = list()
            continue
        if l == "]":
            variant_reading = False
            coded_motifs = new_motifs
            continue
        if variant_reading:
            for m in range(len(coded_motifs)):
                new_motifs.append(coded_motifs[m] + l)
        else:
            for m in range(len(coded_motifs)):
                coded_motifs[m] += l
    return coded_motifs

def read_color_file(color_filename, add_rgb = False):
    if color_filename == None:
        return (dict(), list())
    feature_colors = dict()
    feature_list = list()
    color_file = open(color_filename)
    for string in color_file:
        #{TMHMM}	rgb(255,0,0)
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        fields[0] = fields[0].strip("{}")
        if add_rgb == True:
            fields[1] = "rgb" + fields[1]

        features_to_add = list()
        if fields[0].count("[") != 0: # This is motif and it should be expanded
            features_to_add.extend(vary_motif(fields[0]))
        else:                         # No expansion required
            features_to_add.append(fields[0])
        for f in features_to_add:
            if f in feature_colors:
                print ("Warning: color for '%s' feature was already given: %s" % (f, feature_colors[f]))
                print ("         Changing color to %s" % fields[1])
            feature_colors[f] = fields[1]
            feature_list.append(f)
    color_file.close()
    return (feature_colors, feature_list)

def read_group_file(filename): #FIX: version 1.98 (many columns are allowed now!)
    """
    Assignment of a protein ID to a certain group: should be a file with
    any number of tab-delimited columns, the first one is a protein ID.
    #COG0055_COG1155.ids	Organism_energy_type	c1_specificity_type
    YP_676966.1	                H+	                H+|S		
    ^                           ^                       ^
    Target ID                   Target group 1          Target group 2 <...>
    """
    group = dict()
    group_file = open(filename)
    for string in group_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        if len(fields) > 1:
            seq_id = fields[0]
            match = re.search("([^\.]+)\.\d+$", seq_id)
            if match != None:
                seq_id = match.group(1)
            group[seq_id] = fields[1:]
    group_file.close()
    return group

#------------------------------------------------------------------------------
#                            (2) TMHHM files
#------------------------------------------------------------------------------
def get_feature_from_topology(topology):
    """
    Converts TMHMM topology string to the feature string, i.e.
    returns '[TMHMM] 5..8,10..11' from 'o5-8i10-11o'
    """
    result = "[TMHMM] "
    helix = list()
    fields = topology.strip("oi").split("o")
    for f in fields:
        helix.extend(f.split("i"))
    for h in helix:
        begin = int(h.split("-")[0])
        end = int(h.split("-")[1])
        result += "%i..%i," % (begin, end)
    result = result.strip(",")
    return result

def read_TMHMM_output(input_filename, min_helix_num = 1):
    """
    Reads file <input_filename> (TMHMM output 'one line per protein')
    Returns: dictionary connecting IDs of proteins which have helices to
    corresponding feature strings. 'Having helices' means it has more than
    <min_helix_num> helices.
    """
    TMHMM_result = dict()
    input_file = open(input_filename, "r")
    for string in input_file:
        #sp|P18935|CYB_DROME	len=378	ExpAA=201.00	First60=22.66	PredHel=9	Topology=o34-56i77-99o114-136i141-159o179-201i230-252o289-311i318-340o350-369i
        string = string.strip()
        if len(string) != 0:
            fields = string.split("\t")
            if len(fields) != 6:
                print ("FATAL ERROR: this string in file is not proper TMHMM string:")
                print (string)
                sys.exit()
            curr_id = fields[0]
            helix_num = int(fields[4].split("=")[1])
            if helix_num < min_helix_num:
                continue
            topology = fields[5].split("=")[1]
            TMHMM_result[curr_id] = get_feature_from_topology(topology)
    input_file.close()
    return TMHMM_result

def read_DeepTMHMM_output(input_filename, min_TM_num = 1):
    """
    Reads TMRs.gff3 output of the DeepTMHMM.
    Returns: dictionary connecting IDs of proteins with enought TM elements to
    the corresponding feature strings. 'Enought TM elements' means it has more than
    <min_TM_num> helices or beta-strands.
    """    
    id_to_list_of_TM = dict()
    id_to_list_of_strings = dict()
    id_order = list()
    input_file = open(input_filename)
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] in ["#", "/"]:
            continue
        #1ldf_A	TMhelix	15	32
        #1qd6_C	Beta sheet	12	16
        try:
            fields = string.split("\t")
            protein_id = fields[0]
            element_type = fields[1]
            begin = int(fields[2])
            end = int(fields[3])
            if not protein_id in id_to_list_of_strings:
                id_to_list_of_strings[protein_id] = list()
            id_to_list_of_strings[protein_id].append(string)    
            if element_type in ["TMhelix", "Beta sheet"]:
                if not protein_id in id_to_list_of_TM:
                    id_to_list_of_TM[protein_id] = list()
                    id_order.append(protein_id)
                id_to_list_of_TM[protein_id].append((begin, end))
        except IndexError:
            print("Wrong number of fields in the following line:")
            print(string)
    input_file.close()

    TMHMM_result = dict()
    for protein_id in id_order:
        if len(id_to_list_of_TM[protein_id]) >= min_TM_num:
            id_to_list_of_TM[protein_id].sort(key = lambda k: k[0], reverse = False)
            TMHMM_result[protein_id] = "[TMHMM] "
            for pair in id_to_list_of_TM[protein_id]:
                TMHMM_result[protein_id] += "%i..%i," % (pair[0], pair[1])
            TMHMM_result[protein_id] = TMHMM_result[protein_id].strip(",")
    return (TMHMM_result, id_to_list_of_strings)

#------------------------------------------------------------------------------
#                            (3) Pfam files
#------------------------------------------------------------------------------
class Pfam_domain:
    def __init__(self, name, ac, description):
        self.name = name
        self.ac = ac
        self.description = description

def get_values(field):
    values = field.split("..")
    begin = int(values[0])
    end = int(values[1])
    evalue = float(values[2])
    score = float(values[3])
    try:
        hmm_begin = int(values[4])
        hmm_end = int(values[5])
    except ValueError:
        hmm_begin = 1
        hmm_end = 100500
    return (begin, end, evalue, score, hmm_begin, hmm_end)

def get_feature_from_Pfam(id_to_domains, add_score = False):
     """
     Converts double-hash <id_to_domains> to a single string hash id-to-feature.
     Could add score information (e-value, score, hmm coordinates) to the feature string,
     if <add_score> is True
     """
     Pfam_result = dict()
     for curr_id in id_to_domains.keys():
         feature_string = ""
         for curr_name in id_to_domains[curr_id].keys():
             feature_string += "[%s] " % curr_name
             fields = id_to_domains[curr_id][curr_name].strip().split(" ")
             for f in fields:
                 (begin, end, e_value, score, hmm_begin, hmm_end) = get_values(f)
                 if add_score:
                     feature_string += "%i..%i..%s..%s..%s..%s," % (begin, end, e_value, score, hmm_begin, hmm_end)
                 else:
                     feature_string += "%i..%i," % (begin, end)
             feature_string = feature_string.strip(",")
             feature_string += "\t"
         feature_string.strip("\t")
         Pfam_result[curr_id] = feature_string
     return Pfam_result

def remove_hit(remove_begin, remove_end, remove_name, protein_id, all_domains):
    fields = all_domains[protein_id][remove_name].strip().split(" ")
    all_domains[protein_id][remove_name] = ""
    for f in fields:
        (begin, end, evalue, score, hmm_begin, hmm_end) = get_values(f)
        if (begin == remove_begin) and (end == remove_end): # Removing
            continue
        all_domains[protein_id][remove_name] += "%s..%s..%s..%s..%s..%s " % (begin, end, evalue, score, hmm_begin, hmm_end)
    all_domains[protein_id][remove_name] = all_domains[protein_id][remove_name].strip()

    if all_domains[protein_id][remove_name] == "": # This whole domain is marked for removal
        all_domains[protein_id][remove_name] = "0..0..100..0..0..0"

def check_overlap(target_domain, target_domain_name, protein_id, domains, threshold):
    #print ("  OVERLAP FOR DOMAIN: %s" % target_domain_name)
    target_fields = target_domain.strip().split(" ")
    c = 0
    t = 0
    for tf in target_fields:                               # Itterating domain hits
        (target_begin, target_end, target_evalue, target_score, target_hmm_begin, target_hmm_end) = get_values(tf)
        target_len = target_end - target_begin + 1
        #print ("  TBEGIN = %i, TEND = %i, TEVALUE = %s" % (target_begin, target_end, target_evalue))

        for domain_name in domains[protein_id].keys():
            if domain_name == target_domain_name: # This is the same domain; skipping
                continue
            fields = domains[protein_id][domain_name].strip().split(" ")
            for f in fields:                               # Itterating domain hits
                (curr_begin, curr_end, curr_evalue, curr_score, curr_hmm_begin, curr_hmm_end) = get_values(f)
                #print ("    CBEGIN = %i, CEND = %i, CEVALUE = %s" % (curr_begin, curr_end, curr_evalue))
                curr_len = curr_end - curr_begin + 1
                overlap = min(curr_end, target_end) - max(curr_begin, target_begin) + 1
                if overlap < 0:
                   overlap = 0
                if 100 * float(overlap) / max(target_len, curr_len) > threshold:
                    if curr_evalue < target_evalue:
                        #print ("    Current (%s) better than target (%s) for %s (%s vs. %s)" % (domain_name, target_domain_name, protein_id, curr_evalue, target_evalue))
                        remove_hit(target_begin, target_end, target_domain_name, protein_id, domains)
                        c += 1
                    elif (target_evalue < curr_evalue):
                        #print ("    Target (%s) better than current (%s) for %s (%s vs. %s)" % (target_domain_name, domain_name, protein_id, target_evalue, curr_evalue))
                        remove_hit(curr_begin, curr_end, domain_name, protein_id, domains)
                        t += 1
                    else: # Rare case of e-value equality
                        remove_hit(curr_begin, curr_end, domain_name, protein_id, domains)
                        t += 1
    return (c, t)

def filter_Pfam_hits(id_to_domains, threshold):
    print ("Filtering of domains started...")
    n = 0
    for curr_id in id_to_domains.keys():
        #print ("OVERLAP FOR PROTEIN: %s" % curr_id)
        for curr_name in id_to_domains[curr_id].keys(): # Itterating by domains
            (c, t) = check_overlap(id_to_domains[curr_id][curr_name], curr_name, curr_id, id_to_domains, threshold)
            n += c + t
    print ("Total %i field removements were done!" % n)

    filtered = dict()
    d = 0
    for curr_id in id_to_domains.keys():
        for curr_name in id_to_domains[curr_id].keys():
            fields = id_to_domains[curr_id][curr_name].strip().split(" ")
            remove = False
            for f in fields:
                (begin, end, evalue, score, hmm_begin, hmm_end) = get_values(f)
                if (begin == 0) and (end == 0):
                    d += 1
                    remove = True
            if not remove:
                if not curr_id in filtered:
                    filtered[curr_id] = dict()
                filtered[curr_id][curr_name] = id_to_domains[curr_id][curr_name]
    print ("Total %i domain removements were done!" % d)
    return filtered

def unite_same_Pfam_hits(id_to_domains, max_distance = 50, max_hmm_overlap = 20): #FIX: version 1.96
    """
    Currently the following parameters are used by default:
    <max_distance> = less than 50 aminoacids should be between domain pieces
    <max_hmm_overlap> = less overlap than this portion of the largest HMM parts (20%) between coordinates in HMM
    Change in version 2.9: now scores of domain merges are summed, and e-values are multiplied
    (previously minimal e-value and maximum score of the merged parts were taken).
    Change in version 3.0: regions in HMM should also follow each other, or circular permutation would be not seen!
    """
    n = 0
    for protein_id in id_to_domains.keys():
        #print ("Protein_id = %s" % protein_id)
        for domain_name in id_to_domains[protein_id].keys():
            curr_data = id_to_domains[protein_id][domain_name].strip().split(" ")
            if len(curr_data) > 1: # More than a single hit of this domain found
                curr_data.sort(key = lambda k: int(k.split("..")[0]), reverse = False)
                p = 0
                #print ("  Domain occuring more than once = %s" % domain_name)
                while p < len(curr_data) - 1:
                     (fst_begin, fst_end, fst_evalue, fst_score, fst_hmm_begin, fst_hmm_end) = get_values(curr_data[p])
                     (snd_begin, snd_end, snd_evalue, snd_score, snd_hmm_begin, snd_hmm_end) = get_values(curr_data[p + 1])
                     try:
                         curr_hmm_overlap = min(fst_hmm_end, snd_hmm_end) - max(fst_hmm_begin, snd_hmm_begin) + 1
                         fst_hmm_size = fst_hmm_end - fst_hmm_begin
                         snd_hmm_size = snd_hmm_end - snd_hmm_begin

                         if ((snd_begin - fst_end) < max_distance) and (curr_hmm_overlap < (max(fst_hmm_size, snd_hmm_size) * max_hmm_overlap / 100)):
                             if (fst_hmm_begin < snd_hmm_begin): #FIX: version 3.0 (added condition on HMM coordinates)
                                 #curr_data[p] = "%s..%s..%s..%s..%s..%s" % (fst_begin, snd_end, min(fst_evalue, snd_evalue), max(fst_score, snd_score), fst_hmm_begin, snd_hmm_end)
                                 curr_data[p] = "%s..%s..%s..%s..%s..%s" % (fst_begin, snd_end, fst_evalue * snd_evalue, fst_score + snd_score, fst_hmm_begin, snd_hmm_end) #FIX: version 2.9
                                 curr_data.pop(p + 1)
                                 n += 1
                                 p -= 1
                             #else:
                             #    print ("    Yes! protein_id = %s" % protein_id)
                             #    print (curr_data)
                             #    print (fst_begin, fst_end, fst_evalue, fst_score, fst_hmm_begin, fst_hmm_end)
                             #    print (snd_begin, snd_end, snd_evalue, snd_score, snd_hmm_begin, snd_hmm_end)
                             #    print ("Domain distance: %s; HMM overlap: %s; fst_hmm_size = %s; snd_hmm_size = %s" % (snd_begin - fst_end, curr_hmm_overlap, fst_hmm_size, snd_hmm_size))
                             #    raw_input()
                     except:
                         print ("Bad data found in <unite_same_Pfam_hits> method: '%s'" % curr_data[p])
                         raw_input()
                     p += 1
            id_to_domains[protein_id][domain_name] = " ".join(curr_data)
    print ("Total %i unitings of the same domains were done!" % n)
    return id_to_domains

def read_Pfam_output(input_filename, max_e_value, filter_hits, threshold, add_score = False, unite_same = False, max_distance = 50, max_hmm_overlap = 20, do_not_get_features = False, use_c_evalue = False, hmmsearch_output = False, use_accession = False):
    """
    Method reads <input_filename> produced under -domtblout option by hmmscan. It is
    working with hits with i (independent) e-value not higher than that given <max_e_value>.
    If <filter_hits> is set to True, it will apply filtering on the overlap <threshold>,
    which is currently calculated as:
    100 * float(overlap) / max(target_len, curr_len)
    If <add_score> is True, domain information will be filled with the score.
    If <unite_same> is True, hits of the same domain found (a) near each other (less than <max_distance>
    residues between corresponding regions in a protein) and (b) with different parts of the HMM models
    (less than <max_hmm_overlap> portion of the larger region in HMM should overlap with smaller) will
    be united
    If <do_not_get_features> is True, features are not obtained and <id_to_domains> are returned
    as is
    If <use_c_evalue> if True, conditional e-value (column #11) will be checked against <max_e_value>
    instead of independent e-value (column #12)

    If <hmmsearch_output> is True, query and target fields are changed

    If <use_accession> is True, not domain 'name' (first column) but its 'accession' (second column
    without version mark) will be used.

    Results tuple of two values:
    (1) Hash of protein IDs to a string configured from domains found in a format like this:
    protein_id -> [domain name] s1..e1..evalue..(score)..hmm_s1..hmm_e1,s2..e2..e..(s)..hmm_s2..hmm_e2\t[another domain] <etc> #FIX: version 3.0
    (2) Hash with information about domains found: domain names as keys and <Pfam_domain>
    objects as values
    """
    id_to_domains = dict()  # Hash of hashes: first by protein id, then by domain name.
    domains = dict() # Hash of domain names as keys and <Pfam_domain> objects as values
    input_file = open(input_filename, "r")
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
        #  target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
        #  [0]                    [1]       [2]       [3]               [4]      [5]       [6]   [7]    [8]  [9] [10]  [11]        [12]    [13]  [14]  [15] [16]  [17]  [18]   [19]  [20] [21] [22]
        # ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
        #MatC_N               PF07158.6    149 YP_003148402.1       -            447     1e-56  190.9  19.6   1   2   2.8e-60     1e-56  190.9  13.6     1   149     1   149     1   149 0.99 Dicarboxylate carrier protein MatC N-terminus
        fields = re.split(" +", string, 22)
        try:
            domain_name = fields[0]
            protein_id = fields[3]
            if hmmsearch_output: #FIX: version 3.3
                domain_name = fields[3]
                protein_id = fields[0]
            domain_ac = fields[1].split(".")[0]
            if use_accession:
                domain_name = domain_ac #FIX: version 3.4
                domain_ac = fields[0]
            domain_description = fields[22]
            evalue = None
            if use_c_evalue:
                evalue = fields[11]
            else:
                evalue = fields[12]
            score = fields[13]
            hmm_length = fields[2]
            hmm_begin = fields[15]
            hmm_end = fields[16]
            begin = fields[17]
            end = fields[18]
        except IndexError:
            print ("This string is no standart HMMer domtable output: %i fields detected" % len(fields))
            print (string)
            print (fields)
            sys.exit()
        #print ("NAME = %s, AC = %s, DESCR = %s, PROTEIN = %s, begin = %s, end = %s" % (domain_name, domain_ac, domain_description, protein_id, begin, end))

        #---- 1) E-value check
        if float(evalue) > float(max_e_value):
            continue

        #---- 2) Adding to hash with protein id as keys
        if not protein_id in id_to_domains:
            id_to_domains[protein_id] = dict()
            id_to_domains[protein_id][domain_name] = ""
        else:
            if not domain_name in id_to_domains[protein_id]:
                id_to_domains[protein_id][domain_name] = ""
        hmm_covered = round((int(hmm_end) - int(hmm_begin) + 1) * 100 / int(hmm_length)) # FIX: version 2.1
        id_to_domains[protein_id][domain_name] += "%s..%s..%s..%s..%s..%s..%s " % (begin, end, evalue, score, hmm_begin, hmm_end, hmm_covered)
        #---- 3) Adding to hash with domain names as keys
        if not domain_name in domains:
            curr_domain = Pfam_domain(domain_name, domain_ac, domain_description)
            domains[domain_name] = curr_domain
    input_file.close()
    for pid in id_to_domains.keys(): # FIX: version 3.1 (removement of the last space)
        for domain in id_to_domains[pid].keys():
            id_to_domains[pid][domain] = id_to_domains[pid][domain].strip()
 
    if unite_same == True:
        id_to_domains = unite_same_Pfam_hits(id_to_domains, max_distance, max_hmm_overlap)

    if filter_hits == True:
        id_to_domains = filter_Pfam_hits(id_to_domains, threshold)

    if do_not_get_features:
        Pfam_result = id_to_domains;
    else:
        Pfam_result = get_feature_from_Pfam(id_to_domains, add_score)
    return (Pfam_result, domains)

def get_length(COG_record):
    coordinates = COG_record.split(" ")[1].split(",")
    result = 0
    for c in coordinates:
        result += int(c.split("..")[1]) - int(c.split("..")[0]) + 1
    return result

def read_plain_features(features_filename, no_sort = True):
    """
    Method reads output of the <obtain_plain_features.py> script.
    By default, it would not perform sorting of the results, but
    if <no_sort> is False, sorting would be performed by the total
    length of regions covered by a COG in ascending order (I really
    don't know for what possible reason did I do this)
    """
    id_to_COGs = dict()
    features_file = open(features_filename)
    for string in features_file:
        string = string.strip()
        if len(string) == 0:
            continue
        #UP_000000191	[COG0546] 14..168	[COG1051] 163..286
        fields = string.split("\t")
        protein_id = fields.pop(0)
        if no_sort: #FIX: version 2.9
            fields_sorted = fields
        else:
            fields_sorted = sorted(fields, key=get_length)
        for i in range(len(fields_sorted)):
            fields_sorted[i] = fields_sorted[i].split(" ")[0].strip("[]")
        id_to_COGs[protein_id] = fields_sorted
    features_file.close()
    return id_to_COGs
#------------------------------------------------------------------------------
#                            (4) Other HMMer files
#------------------------------------------------------------------------------
def read_HMMer_table(filename, max_e_value, report_file = None):
    """
    Method reads table output produced under hmmsearch -tblout option from <filename>
    and obtains a dictionary of hit id to hit e-value (if e-value is lower or equals to
    the <max_e_value>).
    If <report_file> is not None, report information will be stored to this file
    (it should be already opened, i.e. this is no filename).
    """
    hits = dict()
    input_file = open(filename)
    for string in input_file:
        string = string.strip()
    #                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
    # target name      accession  query name accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
    #     [0]             [1]        [2]        [3]         [4]     [5]    [6]     [7]     [8]   [9]    [10][11][12][13][14][15][16][17]         [18:]
    # ----------------            ---------- ---------  --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
    #167039788|COG0370  -          database       -      1.9e-283  949.5   9.5  2.3e-283  949.2   6.6   1.0   1   0   0   1   1   1   1 Thermoanaerobacter sp. X514
        if len(string) != 0:
            if string[0] != "#":
                fields = re.split(" +", string)
                if len(fields) < 18:
                    print ("FATAL ERROR: the file %s is not a proper HMMer table ouput; obtained %i fields, expected at least 18" % (filename, len(fields)))
                    print ("String:\n%s" % string)
                    sys.exit()
                hit_id = fields[0]
                hit_evalue_str = fields[4]
                hit_evalue = float(hit_evalue_str)
                hit_score = fields[5]
                hit_description = ""
                if len(fields) >= 18:
                    hit_description = " ".join(fields[18:])
                if hit_evalue <= max_e_value:
                    hits[hit_id] = hit_evalue
                    if report_file != None:
                        report_file.write("%s\t%s\t%s\t%s\n" % (hit_id, "%s*" % hit_evalue_str, hit_score.replace(".", ","), hit_description))
    input_file.close()
    return hits

def read_HMM_report(filename, reverse_order, type_of_id = "GI"):
    """
    Method reads the report file printed by the <read_HMMer_table> method of this package. Returns
    hits as a list of tuples, sorted ascending or descending (if <reverse> is True) by the 'Score'.
    Tuple values order:
    [0] = first part of id (before |): id or gi
    [1] = second part of id (after |): COG
    [2] = e-value
    [3] = Score
    [4] = other part of sequence description (e.g. organism)
    Works also if hits contain no COG data (e.g. have NCBI/URef format like 'gi|...|ref|...').
    In this case COG value will be 'N/A'
    """
    values = list()
    report_file = open(filename)
    for string in report_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        curr_id = fields[0].split("|")[0]
        curr_COG = fields[0].split("|")[1]
        if (fields[0].count("COG") == 0) and (fields[0].count("*") == 0): # No COG data FIX: v.3.1 - either normal or masked)
            if type_of_id == "GI":
                curr_id = fields[0].split("|")[1]
            if type_of_id == "ID":
                curr_id = fields[0].split("|")[3]
            curr_COG = "N/A"
        curr_evalue = fields[1].strip("*")
        curr_score = fields[2].replace(",", ".")
        curr_description = fields[3]
        values.append((curr_id, curr_COG, float(curr_evalue), float(curr_score), curr_description))
    values.sort(key = lambda k: k[3], reverse = reverse_order)
    report_file.close()
    return values

#------------------------------------------------------------------------------
#                            (5) COG methods
#------------------------------------------------------------------------------
class COG_hit:
    def __init__(self, begin, end, COG, protein_length = None):
        self.begin_end = "%s..%s" % (begin, end)
        self.COG_full = COG
        self.protein_length = protein_length

    def get_COG(self):
        return self.COG_full.split("/")[0]

    def set_COG(self, new_COG):
        quality = self.get_problem_code()
        self.COG_full = new_COG
        if quality != None:
            self.COG_full += "/%s" % quality

    def get_problem_code(self):
        if self.COG_full.count("/") != 0:
            return self.COG_full.split("/")[1]
        else:
            return None

    def is_normal(self):
        if self.COG_full.count("/") != 0:
            return False
        else:
            return True

    def get_feature(self):
        return "[%s] %s" % (self.COG_full, self.begin_end)

    def begin(self):
        return self.begin_end.split("..")[0]

    def end(self):
        return self.begin_end.split("..")[1]

    """
    def get_problem_data(self):
        problem_hash = {1 : "Protein > COG", 2 : "Protein < COG", 3 : "Poor alignment"}
        problem_code = self.get_problem.code()
        if problem_code in problem_hash:
            return problem_hash[self.problem_code]
        else:
            return None
    """

def sort_COGs(c1, c2):
    b1 = int(c1.begin_end.split("..")[0])
    b2 = int(c2.begin_end.split("..")[0])

    if b1 < b2:
        return -1
    if b1 > b2:
        return 1
    if b1 == b2:
        return 0

def get_COG_for_gi(COG_assignment, gi):
    result = "N/A"
    if gi in COG_assignment:
        result = ""
        if type(COG_assignment[gi]) == list:
            for cog in COG_assignment[gi]:
                 result += "%s_" % cog.COG_full
            result = result.strip("_")
        else:
            result = COG_assignment[gi].COG_full
    return result

def get_list_of_COGs(COG_assignment, gi):
    """
    Method returns list of COG ids to which given <gi> is attributed in <COG_assignment>
    """
    COGs = list() # List of COGs to which this protein with <gi> ID was attributed
    if gi in COG_assignment:
        #FIX: v.2.7 - now COG hits are sorted according to the position of COG footprint begin!
        COGs_unsorted = list()
        if type(COG_assignment[gi]) != list: # This was not converted to list
            COGs_unsorted.append(COG_assignment[gi])
        else:
            for c in COG_assignment[gi]:
                COGs_unsorted.append(c)
        COGs_sorted = sorted(COGs_unsorted, key = lambda k: k.begin(), reverse = False)
        for c in COGs_sorted:
            COGs.append(c.get_COG())        
    else:
        print ("[WARNING] Method <udav_soft.get_list_of_COGs> did not found protein '%s' in COG assignment!" % gi)
    return COGs 

def write_single_COG_assignment(file_path, COG_assignment, organism):
    output_file = open(file_path, "w")
    for gi in COG_assignment.keys():
        COG_objects = COG_assignment[gi]
        for c in COG_objects:
            curr_string = "%s,%s,%s,%s,%s,%s,%s,\n" % (gi, organism, gi, c.protein_length, c.begin(), c.end(), c.COG_full)
            output_file.write(curr_string)
    output_file.close()

def read_single_COG_assignment(file_path, COG_assignment_large, gi_order, uid_to_org, gi_to_org_file, type_of_id, req_org = None, cog_database_slice_file = None, req_COG = None, slice_only = False):
    """
    Method to read a single file from the COG assignment. It both fills the given
    <COG_assignment_large>, <gi_order> and <uid_to_org> values and returns an assignment
    for this particular file and the last organism name found in the file.

    If <slice_only> is True, data will not be stored in memory - only slice into the required file will be done - FIX: v.3.6
    """
    COG_assignment = dict()
    organism = ""
    single_file = open(file_path)
    #print ("Reading file %s..." % file_path)
    for string in single_file:
        #10957100,Buchnera_aphidicola_APS__Acyrthosiphon_pisum__uid57805,10957100,521,1,521,COG0147/1,
        string = string.strip()
        if len(string) == 0:
            continue
        fields = string.split(",")
        begin = fields[4]
        end = fields[5]
        cog = fields[6]
        protein_length = fields[3]
        if len(fields) == 9: # New format
            #158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,            
            if fields[7] != "0":
                cog = fields[6] + "/" + fields[7]
        if len(fields) == 13: # 2020 format - fragments are merged
            #FLELI_RS14210,GCF_000265505.1,WP_014798679.1,651,1-46=282-651,416,COG1158,COG1158,1,640.5,1.0e-200,422,5-48=49-420
            begin = fields[4].split("-")[0]
            end = fields[4].split("-")[-1]
        gi = None
        if type_of_id == "ID": #FIX: version 1.9
            gi = fields[0]
        elif type_of_id == "locus":
            gi = fields[2]
        #curr_uid = fields[1].split("_")[-1] #FIX: v.2.4 Old-style format of organism name (xxx_uidzzzz) is no longer applicable
        curr_uid = fields[1]
        organism = fields[1]
        if (cog_database_slice_file != None) and (req_org != None):
            if organism in req_org:
                if req_COG != None:
                    if cog.upper() == req_COG.upper():
                        cog_database_slice_file.write("%s\n" % string)  
                else:
                    cog_database_slice_file.write("%s\n" % string)
        if not gi in COG_assignment:
            if gi_to_org_file != None: # FIX: version 1.2
                gi_to_org_file.write("%s\t%s\n" % (gi, fields[1]))
        if not slice_only: #FIX: v.3.6 (extra large COG assignments could be treated)
            if not gi in COG_assignment:
                COG_assignment[gi] = list()
                prot_dict = {"gi" : gi, "org_uid" : curr_uid}
                uid_to_org[curr_uid] = fields[1]
                gi_order.append(prot_dict)
            COG_assignment[gi].append(COG_hit(begin, end, cog, protein_length))
        del fields
    single_file.close()

    for gi in COG_assignment.keys():
        if not gi in COG_assignment_large:
            COG_assignment_large[gi] = list()
        else:
            #This is a weird case!
            print ("WARNING: gi '%s' is found in multiple files!" % gi)
            #sys.exit()
        COG_assignment_large[gi].extend(COG_assignment[gi])
        COG_assignment_large[gi] = sorted(COG_assignment_large[gi], key = lambda cog: cog.begin_end.split("..")[0])
    return (COG_assignment, organism)

def read_COG_assignment(assign_folder, gi_to_org_file = None, type_of_id = "ID", req_org = None, cog_database_slice = None, req_COG = None, large_file = False):
    """
    Method reads contents of the <assign_folder>, which should be an assignment
    of proteins to COGs. If <gi_to_org_file> is given, it will be filled with
    assignment of GI to organism.
    ! <type_of_id> supports upgraded format where two GIs in a string are replaced with
    protein_id (1st) and locus_tag (2nd) FIX: version 1.9
    Returns: 1) dictionary <COG_assignment> with assignment of gi to a list of COG_hit objects;
             2) list of all unique gi as they appear in the file
             3) dictionary of uid to real organism name
    Also supports .cvs file from 2014 release with direct assignment -> FIX: 30.01.2015, version 1.4

    If <req_org> and <cog_database_slice> are not None, a slice of the COG database will be written to
    the file named <cog_database_slice>.

    If <large_file> is True, a slice of the COG database will be made first (without a storage in memory),
    and then it will be red properly FIX: v.3.6
    """
    if (type_of_id != "ID") and (type_of_id != "locus"):
        print ("FATAL ERROR: given type of id '%s' is not supported!" % type_of_id)
        return None
    if (large_file == True) and (cog_database_slice == None):
        print ("FATAL ERROR: large files could only be red with the intermediate slice file!")
        return None
    cog_database_slice_file = None
    if cog_database_slice != None:
        cog_database_slice_file = open(cog_database_slice, "w")
    COG_assignment = dict()
    gi_order = list()
    uid_to_org = dict()
    files = None
    if os.path.isdir(assign_folder):
       files = os.listdir(assign_folder)
    else:
       files = [""]
    n = 0
    for f in files:
        file_path = os.path.join(assign_folder, f).strip("\\/")
        read_single_COG_assignment(file_path, COG_assignment, gi_order, uid_to_org, gi_to_org_file, type_of_id, req_org, cog_database_slice_file, req_COG, large_file)
        n += 1
        print ("Assignment file %i (out of %i) processed!" % (n, len(files)))
    if cog_database_slice != None:
        cog_database_slice_file.close()

    if large_file: #FIX: v.3.6
        read_single_COG_assignment(cog_database_slice, COG_assignment, gi_order, uid_to_org, None, type_of_id, req_org, None, req_COG, False)

    return (COG_assignment, gi_order, uid_to_org)

def read_whog(whog_filename):
    """
    Method for reading assignment of COG to its name. Can work with 'whog' file from the
    initial COG release.
    Also supports .tab file from 2014 release with direct assignment -> FIX: 16.01.2015, version 1.3
    """
    COG_to_description = dict()
    whog_file = open(whog_filename)
    format = "whog"
    split_symbol = " "
    for string in whog_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            format = "tab"
            continue

        if (format == "whog") and (string[0] == "["): # This is a required string in whog file
            #[K] COG1158 Transcription termination factor
            fields = string.split(" ", 2)
            if len(fields) != 3:
                print ("FATAL ERROR: File of the COG-to-name assignment is broken in the string:")
                print (string)
                sys.exit()
            COG_to_description[fields[1]] = (fields[2], fields[0])
        if format == "tab": # This is a direct table (as presented in the 2014 COG release)
           #COG0001	H	Glutamate-1-semialdehyde aminotransferase
           fields = string.split("\t")
           COG_to_description[fields[0]] = (fields[2], "[%s]" % fields[1])
    whog_file.close()
    return COG_to_description

def read_COG_category(category_filename):
    """
    Method reads COG functional category table and returns dictionary like:
    [J] -> Translation, ribosomal structure and biogenesis
    """
    category = dict()
    cat_file = open(category_filename)
    for string in cat_file:
        #J	Translation, ribosomal structure and biogenesis
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t") #FIX: v.2.3 Now can also read 3-column fun-20.tab file (2020 release)
        if len(fields) == 2: # 2014 release
            #J	Translation, ribosomal structure and biogenesis
            category["[%s]" % fields[0]] = fields[1]
        if len(fields) == 3: # 2020 release
            #J	FCCCFC	Translation, ribosomal structure and biogenesis
            category["[%s]" % fields[0]] = fields[2]
        if len(fields) == 4: # 2024 release
            #J	1	FCCCFC	Translation, ribosomal structure and biogenesis
            category["[%s]" % fields[0]] = fields[3]
    cat_file.close()
    return category

def get_COG(assign_folder, cog, output_filename = None, list_filename = None, COG_statistics = None, type_of_id = "ID"):
    """
    Method fully mimics result of the <get_COG.pl> script: prints GIs for the given <cog>
    into the file with <output_filename> (if name given). If <list_filename> is given,
    only GI's belonging to the organisms list will be taken.

    Works with both COG database formats.

    If <COG_statistics> is given, creates a table of all COGs found into the COG assignment,
    including total number of members and number of members in a given <list_filename> (if given).

    Returns: list of all ID (fst, if <type_of_id> is 'ID', or snd, if 'locus> found under
             given criteria (COG and list)
    """
    org_list = dict()
    if list_filename != None:
        org_list = read_plain_list(list_filename)
    statistics_file = None

    gi_dict = dict()
    cog_dict = dict()
    files = os.listdir(assign_folder)
    n = 0
    for f in files:
        file_path = os.path.join (assign_folder, f)
        single_file = open(file_path)
        for string in single_file:
            #10957100,Buchnera_aphidicola_APS__Acyrthosiphon_pisum__uid57805,10957100,521,1,521,COG0147/1,
            string = string.strip()
            if len(string) == 0:
                continue
            fields = string.split(",")
            curr_COG = fields[6].split("/")[0]
            curr_org = fields[1]

            curr_gi = None
            if type_of_id == "ID": #FIX: version 1.9
                curr_gi = fields[0]
            elif type_of_id == "locus":
                curr_gi = fields[2]

            if not curr_COG in cog_dict: # New COG
                cog_dict[curr_COG] = [0, 0]
            cog_dict[curr_COG][0] += 1
            if (list_filename == None) or (curr_org in org_list):
                cog_dict[curr_COG][1] += 1

            if curr_COG == cog:
                if (list_filename == None) or (curr_org in org_list):
                    gi_dict[curr_gi] = True
        single_file.close()

    if output_filename != None:
        output_file = open(output_filename, "w")
        for gi in gi_dict.keys():
            output_file.write("%s\n" % gi)
        output_file.close()

    if COG_statistics != None:
        statistics_file = open(COG_statistics, "w")
        statistics_file.write("#List given: %s\n" % list_filename)
        statistics_file.write("#COG\tOccurence\tIn_list\n")
        for cog in cog_dict.keys():
             statistics_file.write("%s\t%i\t%i\n" % (cog, cog_dict[cog][0], cog_dict[cog][1]))
        statistics_file.close()

    return len(gi_dict.keys())

def check_COG_statistics(input_filename, output_filename, whog_filename, category_filename, min_num, in_list = True):
    """
    Method reads list of COG statistics in <filename> and adds information from files with COG names <whog_filename> and
    functional categories <category_filename>.
    If <in_list> is True, prints information about the COGs which are not represented well in the organisms of the
    list (occurence is less than given <min_num>). Otherwise this information is printed for COGs
    not represented in all release.
    Returns: list of 'good' COGs
    """
    COG_to_description = read_whog(whog_filename)
    category = read_COG_category(category_filename)
    good_COGs = list()
    COG_to_numbers = open(input_filename)
    result = open(output_filename, "w")
    result.write("#G\\B\tCOG\tTotal\tList\tName\tCategory_code\tCategory\n")
    for string in COG_to_numbers:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        cog = fields[0]
        cog_name = "Unk"
        cog_category = "Unk"
        cog_category_full = "Unk"
        try:
            cog_name = COG_to_description[cog][0]
            cog_category = COG_to_description[cog][1]
            if len(cog_category) == 3: #Single category assignment
                cog_category_full = category[cog_category]
            else:
               cog_category_full = ""
               for symbol in cog_category.strip("[]"):
                    curr_category = "[%s]" % symbol
                    cog_category_full += "%s; " % category[curr_category]
               cog_category_full = cog_category_full.strip(" ").strip(";")
        except:
            pass
        total_num = int(fields[1])
        list_num = int(fields[2])
        code = "G"
        if (in_list == True) and (list_num < min_num):
            code = "B"
        if (in_list == False) and (total_num < min_num):
            code = "B"

        if code == "G":
            good_COGs.append(cog)
        result.write("%s\t%s\t%i\t%i\t%s\t%s\t%s\n" % (code, cog, total_num, list_num, cog_name, cog_category, cog_category_full))

    COG_to_numbers.close()
    result.close()
    return good_COGs

def fix_COG_assignment(assign_folder, good_COGs = None, gi_replace = None):
    """
    Method will read assignment of proteins to COGs in <assign_folder> and fix it
    in two possible ways. Be carefull: it will re-write the files, have a backup!
    1) If <good_COGs> is given, it remove all COGs with less members than required:
    only COGs from the <good_COGs> list remain
    2) If a dictionary <gi_replace> is given, it will also replace GIs according to it and
    remove records with 'None' values.
    """
    if (good_COGs == None) and (gi_replace == None):
        print ("COG assignment cannot be fixed, both fixing parametes are 'None'!")
        return
    print ("COG assignment fixing started, %i COGs are considered 'good'" % len(good_COGs))
    files = os.listdir(assign_folder)
    o = 0
    n = 0
    for f in files:
        file_path = os.path.join (assign_folder, f)
        single_file = open(file_path)
        strings = list()
        curr_o = 0
        for string in single_file:
            #10957100,Buchnera_aphidicola_APS__Acyrthosiphon_pisum__uid57805,10957100,521,1,521,COG0147/1,
            string = string.strip()
            if len(string) == 0:
                continue
            fields = string.split(",")
            curr_COG = fields[6].split("/")[0] # Works with both assignment versions

            take_string = True
            if good_COGs != None:         # COG fixing
                if not curr_COG in good_COGs:
                    take_string = False
            if gi_replace != None:        # GI fixing
                if fields[0] in gi_replace:
                    if gi_replace[fields[0]] == None: # deleting
                        take_string = False
                    else:                             # replacing
                        string = string.replace(fields[0], gi_replace[fields[0]])
            if take_string:
                strings.append(string)
            else:
                curr_o += 1
        single_file.close()
        n += 1
        o += curr_o
        print ("Assignment file %i (out of %i) processed; %i strings removed" % (n, len(files), curr_o))
        if curr_o != 0: # Re-writing file
            single_file = open(file_path, "w")
            for string in strings:
                 single_file.write("%s\n" % string)
            single_file.close()
    print ("Total %i strings removed!" % o)

def write_COG_database_per_COG(cog_database, output_dir):
    """
    Method creates multiple .csv files out of given <cog_database>
    corresponding to each separate COG.
    """
    cogdata_filenames = os.listdir(cog_database)
    f = 0
    for cogdata_filename in cogdata_filenames:
        print ("Proceeding file '%s'..." % cogdata_filename)
        cogdata_file = open(os.path.join(cog_database, cogdata_filename))
        for string in cogdata_file:
            #    [0]                                       [1]                                        [2]   [3][4][5]  [6]  [7]
            #NP_603437.1,Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586@GCF_000007325.1_ASM732v1,FN0540,434,1,434,COG0001,0
            string = string.strip()
            if len(string) == 0:
                continue
            if string[0] == "#":
                continue
            fields = string.split(",")
            curr_COG = fields[6]
            COG_filename = os.path.join(output_dir, "%s.csv" % curr_COG)
            if not os.path.isfile(COG_filename): # First occurence of this COG
                f += 1
            COG_file = open(COG_filename, "a")
            COG_file.write("%s\n" % string)
            COG_file.close()
        cogdata_file.close()
    print("DONE! %i files created" % f)

#------------------------------------------------------------------------------
#                            (6) Genbank assembly and Prodigal
#------------------------------------------------------------------------------

def read_assembly_summary(filename):
    """
    This method reads file with the assembly summary provided by NCBI:
    ftp.ncbi.nih.gov:/genomes/genbank/
    """
    return None

def read_protein_table_info(filename, reverse = False):
    """
    This method reads a .table output of the <get_operon.py> script from the
    given <filename>.

    If <reverse> is False, assignment will be as given below. If it is True,
    it would be other way round; protein_id and locus are the keys

    Returns: (1) <GI_to_ID>    - dictionary with GI as keys and protein_id as values
             (2) <GI_to_locus> - dictionary with GI as keys and locus as values
             (3) <gi_dupl>     - list of strings which contain duplicated GIs
             (4) <non_unique>  - list of protein_id and locuses which turned out to be
                                 non-unique
    """
    print ("Assigning GI to ID/locus...")
    GI_to_ID = dict()
    GI_to_locus = dict()
    other_IDs = dict()
    gi_dupl = list()
    non_unique = list()

    print ("    Reading assignment of GI to protein ID / locus...")
    input_file = open(filename, "r") # FIX: 1.95 - no reading of the file into memory
    for string in input_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue

        fields = string.split("\t")
        if len(fields) < 5:
            print ("FATAL ERROR: Not enought elements in string, wrong .table file!")
            print ("Fields: '%s'" % fields)
            sys.exit()
        ID = fields[0]
        gi = fields[1]
        locus = fields[5]
        if gi in GI_to_ID:
            #print ("WARNING: GI '%s' occured multiple times! First assignment (%s) is used" % (gi, other_id))
            gi_dupl.append(string)
            continue
        if ID in other_IDs:
            #print ("FATAL WARNING: ID '%s' is not unique" % ID)
            non_unique.append(ID)
        if locus in other_IDs:
            #print ("FATAL WARNING: locus '%s' is not unique" % ID)
            non_unique.append(locus)
        other_IDs[ID] = True
        other_IDs[locus] = True
        if reverse == True:
            GI_to_ID[ID] = gi
            GI_to_locus[locus] = gi
        else:
            GI_to_ID[gi] = ID
            GI_to_locus[gi] = locus
    input_file.close()
    print ("DONE!")
    return (GI_to_ID, GI_to_locus, gi_dupl, non_unique)

#------------------------------------------------------------------------------
#                            (7) BLAST
#------------------------------------------------------------------------------

class BLAST_hit:
    def __init__ (self, curr_id, start, end, score = None):
        self.curr_id = curr_id
        self.start = int(start)
        self.end = int(end)
        try:
            self.score = float(score)
        except TypeError:
            self.score = None

def read_blast_table(filename, exact_id = True):
    """
    Method returns dictionary of <BLAST_hit> objects with protein_ids of subjects as keys.
    If <exact_id> is True, ids are taken 'as is', and if it is False, they are deduced from URef format
    (e.g., if id was 'gi|Unk|ref|NP_523922.1', only 'NP_523922.1' will be taken)
    """
    result = dict()
    blast_table = open (filename, "r")
    for string in blast_table:
        string = string.strip()
        if len(string) != 0:
            fields = string.split("\t")
            if len(fields) != 12:
                print ("FATAL ERROR: format of blast output is wrong. Expected 12 columns")
                print ("STRING: %s" % string)
                sys.exit()            
            curr_id = fields[1]
            if not exact_id:
                try:
                    parts = curr_id.split("|")
                    curr_id = parts[3]
                except IndexError:
                    pass
            curr_start = fields[8]
            curr_end = fields[9]
            curr_score = fields[11]
            result[curr_id] = BLAST_hit(curr_id, curr_start, curr_end, curr_score)
    blast_table.close()
    return result

class BLAST_hit_detailed:
    def __init__(self, query_id, query_org, query_COG, query_length, target_id, target_org, target_COG, target_length, hit_position, string):
        self.query_id = query_id
        self.query_org = query_org
        self.query_COG = query_COG
        self.query_length = query_length
        self.target_id = target_id
        self.target_org = target_org
        self.target_COG = target_COG
        self.target_length = target_length
        self.hit_position = hit_position
        self.string = string

    def get_data(self):
        fields = self.string.split()
        target_identity = fields[2]
        target_evalue = fields[10]
        alignment_length = fields[3] 
        curr_data = "%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s" % (self.query_id, self.target_id, self.query_org, self.target_org, 
                                                        self.hit_position, target_identity, target_evalue, alignment_length)
        curr_data += "\t%s\t%s" % (self.query_COG, self.target_COG)
        return curr_data

    def get_evalue(self):
        fields = self.string.split()
        target_evalue = float(fields[10])
        return target_evalue

    def get_score(self):
        fields = self.string.split()
        target_score = float(fields[11])
        return target_score

    def get_identity(self):
        fields = self.string.split()
        target_identity = fields[2]
        return target_identity

    def get_query_coverage(self):
        fields = self.string.split()
        result = None
        if self.query_length != 0:
            qstart = int(fields[6])
            qend = int(fields[7])
            result = 100*(qend - qstart + 1)/self.query_length
        return result
        
    def get_target_coverage(self):
        fields = self.string.split()
        result = None
        if self.target_length != 0:
            tstart = int(fields[8])
            tend = int(fields[9])
            result = 100*(tend - tstart + 1)/self.target_length
        return result
