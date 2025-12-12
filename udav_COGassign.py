#!/usr/bin/env python
"""
This module contains basic methods of the <COGassign> program
@ Daria Dibrova aka udavdasha
"""
import sys, os, argparse
import platform
import time
import random
import subprocess
cannot_make_svg = False
try:
    import svgwrite
except ImportError:
    cannot_make_svg = True
    print ("WARNING: Failed to import the module 'svgwrite'")
    print ("Vector scheme will not be generated")

if platform.system() == "Windows":
    sys.path.append("D:\\UdavBackup\\_Complete_genomes\\_scripts")
else:
    sys.path.append("/mnt/Data/scripts")

try:
    import udav_soft
except:
    from . import udav_soft
#========================================================================================
curr_version = 3.0

def read_color_file(filename):
    """
    Required format: <dom_name>\t<color>\n, e.g.:
    COG00055	#ff00ff
    """
    domain_to_color = dict()
    color_file = open(filename)    
    for string in color_file:
        string = string.strip()
        if len(string) == 0:
            continue
        domain = string.split("\t")[0]
        color = string.split("\t")[1].lower()
        domain_to_color[domain] = color
    color_file.close()
    return domain_to_color

def get_random_color():
    r = random.randint(128, 255)
    g = random.randint(128, 255)
    b = random.randint(128, 255)
    random_color = "#%02x%02x%02x" % (r, g, b)
    return random_color

def get_different_colors():
    colors = list()
    red   = [0,   0,   255,   0, 255, 255]
    green = [0,   255, 0,   255, 0,   255]
    blue  = [255, 0,   0,   255, 255,   0]   
    for i in range(6):
        r = red[i]
        g = green[i]
        b = blue[i]    
        colors.append("#%02x%02x%02x" % (r, g, b))
    for i in range(6):
        r = round(red[i] / 2)
        g = round(green[i] / 2)
        b = round(blue[i] / 2)
        colors.append("#%02x%02x%02x" % (r, g, b))  
    for i in range(6):
        r = round(red[i] / 4)
        g = round(green[i] / 4)
        b = round(blue[i] / 4)
        colors.append("#%02x%02x%02x" % (r, g, b))        
    for i in range(6):
        r = 64 + round(red[i] / 2)
        g = 64 + round(green[i] / 2)
        b = 64 + round(blue[i] / 2) 
        colors.append("#%02x%02x%02x" % (r, g, b))        
    return colors

def get_proteins(strings, unique_required = False, quiet = False):
    """
    Method reads proteins in 'Basic' format and also checks if identical
    protein IDs were detected (this would result in erroneous domain assignment)
    if <unique_required> option is True. If <quiet> is True, it will not
    print anything to the console.
    """
    proteins = list()
    non_unique = None
    for string in strings:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == ">":
            protein_id = string.split(" ", 1)[0].strip(">")
            proteins.append([protein_id, "", list(), string])
        else:
            proteins[-1][1] += string
    if unique_required:
        u = 0
        n = 0
        ids_detected = dict()
        non_unique = dict()
        if not quiet: print ("Non-unique proteins will be filtered")
        while u < len(proteins):
            protein_id = proteins[u][0]
            if protein_id in ids_detected: # Was already detected
                if not quiet: print ("[WARNING]: protein ID '%s' appears more than once, only the first occurence counts!" % protein_id)
                non_unique[protein_id] = True
                n += 1
                proteins.pop(u)
                u -= 1
            else:
                ids_detected[protein_id] = True
            u += 1
        if not quiet: print ("Total %i non-unique protein IDs removed; number of proteins remained: %i" % (n, len(proteins)))
    return (proteins, non_unique)

def get_filtered_domains(Pfam_result, all_domains):
    filtered_domains = dict()
    for protein_id in Pfam_result.keys():
        string = "%s\t%s" % (protein_id, Pfam_result[protein_id])
        domain_dict = get_domains(string)
        for domain in domain_dict.keys():
            filtered_domains[domain] = all_domains[domain]
    return filtered_domains

def get_domains(string):
    """
    Method obtains a single string from the features file (produced by <obtain_features.py>
    script) and returns <domains> dictionary for this protein: keys are domain names,
    values are coordinates
    """
    fields = string.strip().split("\t")
    protein_id = fields.pop(0)
    domains = dict()     
    for f in fields:
        name = f.split(" ")[0].strip("[]")
        coordinates = f.split(" ")[1]
        domains[name] = coordinates
    return domains

def sort_domains(domains):
    """
    Method sorts domains obtained via <get_domains> method by their length in descending order,
    i.e. the most long will be located at the top.
    Returns a list of tuples: (domain_name, region). If the same domain name occures multiple times, it will occur
    multiple times in this list
    """
    domain_list = list()
    for name in domains.keys():            
        regions = domains[name].split(",") # If this domain is found in protein multiple times
        for region in regions:
            domain_list.append((name, region))
    domain_list.sort(key = lambda k: int(k[1].split("..")[1]) - int(k[1].split("..")[0]), reverse = True)
    return domain_list

def run_hmmscan(input_file, out_prefix, hmm_database, path_to_hmmscan = None):
    """
    Method runs hmmscan with given parameters and returns a name of resulting file
    (<out_prefix>.table.domains) and <full_proteins> list of all proteins in the
    input FASTA file.
    <input_file>         - it is expected to contain protein sequences in FASTA
                           format (could be a real file object or a <sys.stdin>)
    <out_prefix>         - a prefix for current output files:
                           <out_prefix>.domains - main output with domain
                                                  coordinates and scores
                           <out_prefix>.legend  - legend for coloring and
                                                  domain description
    <hmm_database>       - path to the hmm database file (e.g., COG_database.hmm
                           or Pfam-A.hmm)
    """
    curr_entry = "".join(input_file.readlines())
    if curr_entry.count(">") == 0:
        raise IOError("[..FATAL ERROR..] Input does not contain fasta-formatted sequences!")
    strings = curr_entry.split("\n")
    (proteins, non_unique) = get_proteins(strings, True)
    #------------------------------------------ 1) Writing temp fasta file        
    temp_filename = "%s.fasta" % out_prefix
    temp_file = open(temp_filename, "w")
    for p in proteins: #FIX: version 1.6 (filtering of the non-unique proteins ON)
        temp_file.write(">%s\n" % p[0])
        temp_file.write("%s\n\n" % p[1])
    temp_file.close()
    #------------------------------------------ 2) Running hmmscan
    hmmscan_name = "hmmscan"
    if platform.system() == "Windows":           
        hmmscan_name += ".exe"
    if path_to_hmmscan != None:
        hmmscan_name = os.path.join(path_to_hmmscan, hmmscan_name)

    if not os.path.isfile(hmmscan_name):
        raise IOError("[..FATAL ERROR..] hmmscan was not detected in '%s'" % hmmscan_name)
    result_table = "%s.table.domains" % out_prefix
    result_txt = "%s.txt" % out_prefix
    if os.path.isfile(result_table):
        print ("\nhmmscan have already been run: '%s'" % result_table)
    else:
        print ("\nAnalysis of the sequence started; running <hmmscan>...")
        try:
            hmmscan_output = subprocess.check_output([os.path.join(os.getcwd(), hmmscan_name), "--domtblout", result_table, "-o", result_txt, hmm_database, temp_filename], stderr=subprocess.STDOUT)
            os.remove(result_txt)
        except subprocess.CalledProcessError as error:
            error_filename = "%s.error" % out_prefix
            error_file = open(error_filename, "w")
            error_file.write(str(error.output.decode("utf-8")))
            error_file.close()
            os.remove(temp_filename)
            os.remove(result_table)
            os.remove(result_txt)
            raise IOError ("[..FATAL ERROR..] hmmscan failed to terminate properly! See '%s' file for the error log" % error_filename)
    #------------------------------------------ 3) Temp file removal
    os.remove(temp_filename)
    (full_proteins, non_unique) = get_proteins(strings, False)
    return (result_table, full_proteins)

def main(result_table, filter_thresh, evalue_thresh, unite_same_domains, max_distance, max_hmm_overlap, use_c_evalue = False):
    """
    <result_table>       - name of file with the domtblout result of hmmscan

    <filter_thresh>      - if an overlap between two domains D1 and D2 exceeds
                           this %% of the largest length, i.e.
                           max (len(D1), len(D2)), the one with worse e-value
                           will be removed (filtred). For example, make 0 to filter
                           any overlapping domains (NOT RECOMMENDED, because a
                           small overlap on the tips of domains is OK)

    <evalue_thresh>      - only hits with independent (i) e-value less than this
                           threshold are considered for showing (if an option
                           <use_c_evalue> is True, conditional (c) e-value will be
                           used)

    <hmm_database>       - path to the hmm database file (e.g., COG_database.hmm
                           or Pfam-A.hmm)

    <unite_same_domains> - if True, parts of domain found by different regions of the HMM model
                           should be united if they are found (a) near each other (less than <max_distance>
                           residues between corresponding regions in a protein) and (b) with 
                           different parts of the HMM models (less than <max_hmm_overlap> portion of 
                           the larger region in HMM should overlap with smaller)
    """
    #------------------------------------------ 4) Converting HMMscan results into a proper format

    #sys.stdout = open(os.devnull, "w") # Making STDOUT to be null-device (empty)
    (Pfam_result, domains) = udav_soft.read_Pfam_output(result_table, evalue_thresh, True, filter_thresh, add_score=True, unite_same=unite_same_domains, max_distance = max_distance, max_hmm_overlap = max_hmm_overlap, do_not_get_features=False, use_c_evalue=use_c_evalue) # Pfam result is filtered and e-value/score is added; the same domains could be filtered
    #sys.stdout = sys.__stdout__ # Changing STDOUT back to normal
    #os.remove(result_table)
    domains = get_filtered_domains(Pfam_result, domains)
    return (Pfam_result, domains)

def obtain_data(myargs, proteins, Pfam_result, domains, output_file, legend_file, domain_to_color_preset = None, reqdomains = None, not_strict_req = True):
    #---- 1) Counting domain number and filling proteins[i][2] dictionary
    domain_count = dict()
    for i in range(len(proteins)):
        curr_id = proteins[i][0]
        if curr_id in Pfam_result:
            string = "%s\t%s" % (curr_id, Pfam_result[curr_id])
            domain_dict = get_domains(string)
            proteins[i][2] = sort_domains(domain_dict)                
            for domain in domain_dict.keys():
                if not domain in domain_count:
                    domain_count[domain] = 0
                domain_count[domain] += 1
    #---- 2) Legend printing
    legend_file.write("#Legend created at: %s\n" % time.ctime())
    legend_file.write("#Number\tDomain ID\tColor\tDescription\n")
    color_dict = dict()
    free_color_list = get_different_colors()
    domains_sorted_by_occur = list(domain_count.keys())
    domains_sorted_by_occur.sort(key = lambda k: domain_count[k], reverse = True)
    domain_to_color = dict()
    legend_strings = list()
    for domain in domains_sorted_by_occur:
        curr_color = None 
        if domain_to_color_preset != None: # File with colors for certain domains was given       
            if domain in domain_to_color_preset:
                curr_color = domain_to_color_preset[domain]
                for i in range(len(free_color_list)):
                    if free_color_list[i] == curr_color:
                        free_color_list.pop(i)
                        break
        if curr_color == None:        
            if len(free_color_list) != 0: # There are still available colors
                curr_color = free_color_list[0]
                free_color_list.pop(0)
            else:
                curr_color = get_random_color()
        color_dict[domain] = curr_color
        legend_string = "%s\t%s\t%s\t%s" % (domain_count[domain], domain, curr_color, domains[domain].description)
        legend_strings.append(legend_string)
        legend_file.write("%s\n" % legend_string)
        domain_to_color[domain] = curr_color
    #---- 3) Main output printing
    evalue_used = "i-evalue"
    if myargs["c_evalue"] == True:
        evalue_used = "c-evalue"
    output_file.write("#Output file of <udav_COGassign.py> script created at: %s\n" % time.ctime())
    output_file.write("#Parameters used for the domain search filtering:\n")
    output_file.write("#                         E-value threshold: %s\n" % myargs["evalue_thresh"])
    output_file.write("#                              E-value used: %s\n" % evalue_used)
    output_file.write("#                         Overlap threshold: %.0f%%\n" % myargs["filter_thresh"])
    output_file.write("#                          Domains required: %s\n" % myargs["reqdomains"])
    output_file.write("#                    Unite partial HMM hits: %s\n" % myargs["unite"])
    if myargs["unite"] == True:
        output_file.write("#     Max distance between parts in protein: %s\n" % myargs["max_dist"])
        output_file.write("#     Max partion of larger HMM hit overlap: %s\n" % myargs["max_hmm_overlap"])
    output_file.write("#\n")
    output_file.write("#Protein ID\tLength(aa)\tDomain data\tDomain formula\n")
    for p in proteins:
        #p[2]: [('COG1155', '4..580..5.3e-265..878.3'), ('COG1156', '7..103..6.4e-11..39.5')]
        required_domains = dict()
        if reqdomains != None: #FIX: v.2.0
            for reqdomain in reqdomains.split(" "):
                 required_domains[reqdomain] = True
        domain_data = ""
        for domain in p[2]:
            domain_data += "%s|%s " % (domain[0], domain[1])
            if not_strict_req: # Domain name must not be exactly the same as one in the list, but could contain part or it (e.g. 'Kelch_1' will match 'Kelch')
                for d in required_domains:
                    if d in domain[0]:
                        required_domains.pop(d)
                        break
            else:
                if domain[0] in required_domains:
                    required_domains.pop(domain[0])
        if (reqdomains != None) and (len(required_domains.keys()) != 0): # Not all required domains were found
            continue
        domain_data = domain_data.strip()
        domain_formula = "" #FIX: v.1.8
        sorted_domain_data = sorted(p[2], key = lambda k: int(k[1].split("..")[0]))
        for domain in sorted_domain_data:
            domain_formula += "%s " % (domain[0])
        domain_formula = domain_formula.strip() 
        output_file.write("%s\t%s\t%s\t%s\n" % (p[0], len(p[1]), domain_data, domain_formula))
    output_file.write("# [ok]\n")
    return (domain_to_color, legend_strings)            

def write_svg_scheme(output_filename, proteins, domain_to_color, legend_strings, domain = None, natural = False, add_measure = None):
    """
    Starting with version 1.9, option <natural> can be used to show natural position of domains
    on a protein sequence, like showing them with correct lengths and inter-domain distances.

    <proteins> is the list of lists, each representing data for single protein sequence.
    proteins[i]: [protein_id, protein_sequence, domain_data, fasta_string]
                 [0] protein_id      : part of fasta-formatted sequence name before the first space
                 [1] protein_sequence: complete sequence of this protein
                 [2] domain_data     : list of domains sorted by domain length, e.g. 
                     [('COG1155', '4..580..5.3e-265..878.3'), ('COG1156', '7..103..6.4e-11..39.5')]
                 [3] fasta_string     : fasta-formatted sequence name

    Starting from version 3.0, returns an svg string and if <output_filename> is None, omitts
    writing file.
    """
    if add_measure != None:
        if not add_measure in ["score", "evalue"]:
            raise ValueError("Type of measure to add to the svg scheme '%s' is wrong, should be either 'score' or 'evalue'!" % add_measure)
    # 0) Setting default parameters
    print ("Selected domain to align the picture to it: '%s'" % domain)
    text_style = "font-size:%ipx; font-family:%s" % (22, "Courier New")
    measure_text_style = "font-size:%ipx; font-family:%s font-weight:bold" % (11, "Courier New")
    letter_w = 7.195 * 2 #7.2
    letter_h = 9.6 * 2
    domain_size_x = letter_w * 3
    domain_size_y = letter_h
    between_lines_spacer = 3
    spacer = 2
    field = 5
    max_x_size = 0
    max_y_size = 0
    #max_score=float(800) #fill-opacity:1
    # 1) Printing domains
    #my_svg = svgwrite.Drawing(filename = output_filename, size = ("600px", "1000px"))
    my_svg_main = svgwrite.Drawing(filename = output_filename, size = ("100%", "100%"))
    my_svg = my_svg_main.g()
    min_x = 0 # This value will be used to transform canvas
    for i in range(len(proteins)):
        #proteins[i][2]: [('COG1155', '4..580..5.3e-265..878.3'), ('COG1156', '7..103..6.4e-11..39.5')]
        domain_number = len(proteins[i][2])
        dom_fix = 0
        start_dom_fix = 0
        proteins[i][2].sort(key = lambda k: int(k[1].split("..")[0]), reverse = False)
        if domain != None: # 'Alignment' is required
            for d in range(domain_number):
                if proteins[i][2][d][0] == domain:
                    dom_fix = d
                    start_dom_fix = int(proteins[i][2][d][1].split("..")[0]) # start coordinate of first occurence of this domain
                    break
        
        if natural == True: # Protein length and domain size should be natural FIX: 1.9
            line_length = len(proteins[i][1])
            x = field - start_dom_fix
        else:
            line_length = (domain_number * domain_size_x) + ((domain_number + 1) * spacer)
            x = field - (dom_fix * (domain_size_x + spacer))
        min_x = min(x, min_x) # This value will be used to transform canvas

        line_color = "black"
        y = field + (i * (letter_h + between_lines_spacer))
        protein_line = my_svg_main.line(start = (x, y + letter_h/2),
                              end = (x + line_length, y + letter_h/2),
                              stroke = line_color,
                              fill = "none")
        my_svg.add(protein_line)
        if x + line_length + field > max_x_size:
            max_x_size = x + line_length + field
        if y + letter_h + field > max_y_size:
            max_y_size = y + letter_h + field

        for j in range(domain_number):
            curr_score = float(proteins[i][2][j][1].split("..")[3])
            curr_evalue = proteins[i][2][j][1].split("..")[2]

            score_pos = 1
            #score_pos = curr_score/max_score
            #if score_pos > 1: score_pos = 1
            #if score_pos < 0.1: score_pos = 0.1

            curr_domain = proteins[i][2][j][0]
            curr_color = domain_to_color[curr_domain]
            if natural == True: # Protein length and domain size should be natural FIX: 1.9
                this_domain_start = int(proteins[i][2][j][1].split("..")[0])
                this_domain_end = int(proteins[i][2][j][1].split("..")[1])
                domain_x = x + this_domain_start
                domain_size_x = this_domain_end - this_domain_start + 1
            else:
                domain_x = x + spacer + (j * (domain_size_x + spacer))
            domain_rect = my_svg_main.rect(insert = (domain_x, y),
                              size = ("%ipx" % domain_size_x, "%ipx" % domain_size_y),
                              rx = 10, ry = 10,
                              stroke = "none",
                              fill = curr_color,
                              opacity = score_pos)
            my_svg.add(domain_rect)
            if add_measure != None:
                measure = None
                if add_measure == "score":
                    measure = curr_score
                else:
                    measure = curr_evalue
                measure_text = my_svg_main.text(measure, insert=(domain_x + (domain_size_x/2), y + domain_size_y - 2 - 4), fill="white", style=measure_text_style)
                my_svg.add(measure_text)
    # 2) Writing names of proteins
    longest_name = 0
    for i in range(len(proteins)):
        protein_name = proteins[i][0]
        longest_name = max(len(protein_name), longest_name)
        x = max_x_size
        y = field + (i * (letter_h + between_lines_spacer))
        name_text = my_svg_main.text(protein_name, insert=(x, y + domain_size_y - 2), fill="black", style=text_style)
        my_svg.add(name_text)
    hugest_x = max_x_size + (longest_name * letter_w)
    # 3) Writing legend
    max_y_size = max_y_size + 20
    max_num_size = 3
    max_name_size = 0
    hugest_y = max_y_size
    for s in range(len(legend_strings)):
        curr_name = legend_strings[s].split("\t")[1]
        if len(curr_name) > max_name_size:
            max_name_size = len(curr_name)
    legend_style = "font-size:%ipx; font-family:%s" % (18, "Courier New")
    legend_letter_w = letter_w / 1.22
    for s in range(len(legend_strings)):
        #4	COG0056	#0000ff	FoF1-type ATP synthase, alpha subunit
        fields = legend_strings[s].split("\t")
        curr_color = fields.pop(2)
        text = "\t".join(fields)
        x = field
        y = max_y_size + (s * (letter_h + between_lines_spacer))
        hugest_y = y
        legend_rect = my_svg_main.rect(insert = (x, y),
                          size = ("%ipx" % letter_h, "%ipx" % letter_h),
                          stroke = "none",
                          fill = curr_color)
        my_svg.add(legend_rect)
        num_text = my_svg_main.text(fields[0], insert=(x + field + letter_h, y + letter_h - 2), fill="black", style=legend_style)
        domain_text = my_svg_main.text(fields[1], insert=(x + field + letter_h + (max_num_size * letter_w), y + letter_h - 2), fill="black", style="%s; font-weight:bold" % legend_style)
        descr_text_x = x + field + letter_h + (max_num_size + max_name_size) * letter_w
        descr_text = my_svg_main.text(fields[2], insert=(descr_text_x, y + letter_h - 2), fill="black", style=legend_style)
        hugest_x = max(hugest_x, descr_text_x + (len(fields[2]) * legend_letter_w))
        my_svg.add(num_text)
        my_svg.add(domain_text)
        my_svg.add(descr_text)

    my_svg["transform"] = "translate(%i,0)" % (0 - min_x + field)
    #for element in my_svg.elements:
    #    print (element.__dict__)
    my_svg_main.__setitem__("height", "%spx" % (hugest_y + letter_h + between_lines_spacer))
    my_svg_main.__setitem__("width", "%spx" % (0 - min_x + hugest_x + 5))
    my_svg_main.add(my_svg)
    if output_filename != None:
        my_svg_main.save(pretty = True, indent = 4)
    return my_svg_main.tostring()

def check_and_transform_argument(argument, default, method_to_transform, error_message, bottom_threshold = None, top_threshold = None):
    transformed = default
    if argument != None:
       try:
           transformed = method_to_transform(argument)
           if bottom_threshold != None:
               transformed = max(bottom_threshold, transformed)
           if top_threshold != None:
               transformed = min(top_threshold, transformed)
       except ValueError:
           raise ValueError(error_message)
    return transformed

if __name__ == '__main__':    
    print ("Module <udav_COGassign> (version %s) was executed as a script" % curr_version)
    print ("It will assign proteins from given file to COGs")
    parser = argparse.ArgumentParser(description = "This script will assign proteins from given file to COGs. \
                                                    Current version is %s" % curr_version)
    parser.add_argument("-i", help = "Name of the file with input sequences (FASTA format)", required = False, dest = "input_file")
    parser.add_argument("-o", help = "Name of the main output file", required = True, dest = "output_file")
    parser.add_argument("-l", help = "Name of the legend output file", required = True, dest = "legend_file")
    parser.add_argument("-p", help = "Path to the HMM database (Pfam-A.hmm or other)", required = True, dest = "hmm_database")
    parser.add_argument("-e", help = "E-value threshold for domain filtering", required = False, dest = "evalue_thresh")
    parser.add_argument("-f", help = "Overlap for filtering procedure (%)", required = False, dest = "filter_thresh")
    parser.add_argument("--unite", help = "Unite same domains following each other (use carefully)", action = "store_true", dest = "unite")
    parser.add_argument("--max_dist", help = "If same domains are united, this is a maximum distance between them in a protein", required = False, dest = "max_dist")
    parser.add_argument("--max_hmm_overlap", help = "If same domains are united, this is a maximal overlap between them in HMM coordinates (in % of length of longer HMM hit)", required = False, dest = "max_hmm_overlap")
    parser.add_argument("-c", help = "Filename with pre-defined colors for certain domains", required = False, dest = "color_filename")
    parser.add_argument("-d", help = "The domain which should be 'aligned'", required = False, dest = "domain")
    parser.add_argument("--nosvg", help = "Use this to skip writing both SVG outputs", action = "store_true", dest = "nosvg")
    parser.add_argument("--reqdomains", help = "Use this to print only data for proteins with specified domains, e.g. 'BTB BACK Kelch'", required = False, dest = "reqdomains")
    parser.add_argument("--cevalue", help = "Use this to use conditional e-value (c-evalue) instead of the independent e-value (i-evalue) which is used by default", action = "store_true", dest = "c_evalue")
    parser.add_argument("--hmmscan", help = "Path to the directory with the hmmscan binary file", required = False, dest = "path_to_hmmscan")


    myargs = parser.parse_args()
    #----====---- 1) Arguments check    
    myargs.filter_thresh = check_and_transform_argument(myargs.filter_thresh, 5, int, "FATAL ERROR: Non-int value under -f option '%s'" % myargs.filter_thresh)
    myargs.evalue_thresh = check_and_transform_argument(myargs.evalue_thresh, "1e-05", float, "FATAL ERROR: Non-float value under -e option '%s'" % myargs.evalue_thresh)
    myargs.max_dist = check_and_transform_argument(myargs.max_dist, 50, int, "FATAL ERROR: Non-int value under --max_dist option '%s'" % myargs.max_dist)
    myargs.max_hmm_overlap = check_and_transform_argument(myargs.max_hmm_overlap, 30, int, "FATAL ERROR: Non-int value under --max_hmm_overlap option '%s'" % myargs.max_hmm_overlap)

    if not os.path.isfile(myargs.hmm_database):
        print ("FATAL ERROR: HMM database '%s' does not exist!" % myargs.hmm_database)
        sys.exit()
    #----====---- 2) Checking input and output flows
    input_file = sys.stdin
    output_file = sys.stdout
    out_prefix = "temp"
    close_input = False

    if myargs.input_file != None: # Filename was given
        try:
            real_file = open(myargs.input_file, "r")
            input_file = real_file
            close_input = True
        except:
            print ("FATAL ERROR: file '%s' does not exist!" % myargs.input_file)
            sys.exit()
    else:
        if sys.stdin.isatty():
            print ("FATAL ERROR: <STDIN> is empty, and no input file provided!")
            sys.exit()
    
    output_file = open(myargs.output_file, "w") 
    legend_file = open(myargs.legend_file, "w") 
    out_prefix = myargs.output_file

    t_start = time.time()
    (result_table, proteins) = run_hmmscan(input_file, out_prefix, myargs.hmm_database, myargs.path_to_hmmscan)
    (Pfam_result, domains) = main(result_table, myargs.filter_thresh, myargs.evalue_thresh, myargs.unite, myargs.max_dist, myargs.max_hmm_overlap, myargs.c_evalue)
    domain_to_color_preset = None
    if myargs.color_filename != None:
        domain_to_color_preset = read_color_file(myargs.color_filename)
    (domain_to_color, legend_strings) = obtain_data(myargs.__dict__, proteins, Pfam_result, domains, output_file, legend_file, domain_to_color_preset, myargs.reqdomains)
    t_end = time.time()
    print ("[..SUCCESS..] DONE! Domain search for %i proteins was done in %.1f s" % (len(proteins), (t_end - t_start)))

    if (not myargs.nosvg) and (not cannot_make_svg):
        write_svg_scheme(out_prefix + ".svg", proteins, domain_to_color, legend_strings, myargs.domain)
        write_svg_scheme(out_prefix + ".natural.svg", proteins, domain_to_color, legend_strings, myargs.domain, True)

    if close_input: input_file.close()
    output_file.close()
    legend_file.close()