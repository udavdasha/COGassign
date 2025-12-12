#!/usr/bin/env python
import argparse
try:
    import udav_COGassign_utils
except ImportError:
    from . import udav_COGassign_utils
try:
    import udav_COGassign
except ImportError:
    from . import udav_COGassign
#========================================================================================
curr_version = 1.0
parser = argparse.ArgumentParser(description =
"This script enables an interface between <udav_COGassign.py> and its web server version. \
Current version is %s" % curr_version
)
parser.add_argument("--work_dir", help = "Name of the main working directory", required = True, dest = "work_dir")
parser.add_argument("--jobid", help = "Name of current job id", required = True, dest = "jobid")
parser.add_argument("--overlap", help = "Overlap for filtering procedure (%)", required = True, dest = "overlap")
parser.add_argument("--evalue", help = "E-value threshold for domain filtering", required = True, dest = "evalue")
parser.add_argument("--unite", help = "Unite same domains following each other (use carefully)", action = "store_true", dest = "unite")
parser.add_argument("--max_distance", help = "If same domains are united, this is a maximum distance between them in a protein", required = False, default = "50", dest = "max_distance")
parser.add_argument("--max_hmm_overlap", help = "If same domains are united, this is a maximal overlap between them in HMM coordinates (in % of length of longer HMM hit)", required = False, default = "30", dest = "max_hmm_overlap")
parser.add_argument("--evalue_type", help = "Use 'c' if conditional e-value (c-evalue) should be used or 'i' to use an independent e-value (i-evalue)", required = True, dest = "evalue_type")

myargs = parser.parse_args()
#========================================================================================

def run_udav_COGassign_methods(work_dir, jobid, filter_thresh, evalue_thresh, unite, max_dist, max_hmm_overlap, evalue_type):
    """
    Method is used to produce a file *.domains + *.legend. See description of
    the <udav_COGassign.main> method for other parameters specifications.
    Produces a list of 'proteins' (each 'protein' is a list of 4 values: protein_id, sequence, a list with
    the domain data and a full string with this protein description).
    Returns a name of the final output file (*.domains).
    """
    output_filename = udav_COGassign_utils.get_filename(work_dir, jobid, "domains")
    output_file = open(output_filename, "w")
    fasta_filename = udav_COGassign_utils.get_filename(work_dir, jobid, "fasta")
    domtblout_filename = udav_COGassign_utils.get_filename(work_dir, jobid, "domtblout")
    strings = udav_COGassign_utils.get_strings_from_fasta_file(fasta_filename)
    (proteins, non_unique) = udav_COGassign.get_proteins(strings, True, True)
    if evalue_type == "i":
        c_evalue = False
    else:
        c_evalue = True
    (Pfam_result, domains) = udav_COGassign.main(domtblout_filename, filter_thresh, evalue_thresh, unite, max_dist, max_hmm_overlap, c_evalue)
    args = {"evalue_thresh": evalue_thresh, "c_evalue": c_evalue, "filter_thresh": filter_thresh, "reqdomains" : None, "unite": unite, "max_dist": max_dist, "max_hmm_overlap":max_hmm_overlap}
    (domain_to_color, legend_strings) = udav_COGassign.obtain_data(args, proteins, Pfam_result, domains, output_file, output_file)
    output_file.close()

run_udav_COGassign_methods(myargs.work_dir, myargs.jobid, int(myargs.overlap), float(myargs.evalue), myargs.unite, int(myargs.max_distance), int(myargs.max_hmm_overlap), myargs.evalue_type)
