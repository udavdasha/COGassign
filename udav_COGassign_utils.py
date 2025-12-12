"""
This module contains interface methods for the web-server version of Domain Analyser
(udav_COGassign.py).
@ Daria Dibrova aka udavdasha
"""
curr_version = "1.0"
import sys, os, re
import json
import subprocess
import codecs
try:
    import udav_COGassign
except ImportError:
    from . import udav_COGassign

def get_filename(work_dir, jobid, file_type, filename_only = False):
    """
    This method returns a proper full path to the file of required
    <file_type>. This must be one of the following, otherwise None will be returned:
     * 'fasta' - file with the input protein sequences (under unique IDs only)
     * 'domtblout' - --domtblout output of the hmmscan
     * 'out'- basic output (-o) of the hmmscan
     * 'hmmscan_stderr' - error log of the hmmscan
     * 'hmmscan_stdout' - output log of the hmmscan (should be empty)
     * 'hmmscan_command' - file with the string used to run hmmscan
     * 'domains' - resulting file of the <udav_COGassign.py> (with filtered domains and legend strings)
     * 'svg' - svg-file produced with the <udav_COGassign.py>
     * 'params_domalyser' - file with current parameters for the domain analyser
     * 'database_name' - file with the name of current database
    If <filename_only> is True, no absolute path but a filename only is returned
    """
    result = None
    filenames = {"fasta" : "%s.fasta" % jobid,
                 "domtblout": "%s.domtblout.txt" % jobid,
                 "out": "%s.out.txt" % jobid,
                 "hmmscan_stderr": "%s.hmmscan.stderr.txt" % jobid,
                 "hmmscan_stdout": "%s.hmmscan.stdout.txt" % jobid,
                 "hmmscan_command": "%s.hmmscan.command.txt" % jobid,
                 "domains": "%s.domains" % jobid,
                 "svg"    : "%s.svg" % jobid,
                 "params_domalyser" : "%s.params_domalyser" % jobid,
                 "database_name" : "%s.database_name" % jobid,
                 "id_order"      : "%s.id_order" % jobid,
                 "domalyser_command": "%s.domalyser.command.txt" % jobid,
                 "non_unique" : "%s.non_unique" % jobid
                }
    if file_type in filenames:
        if filename_only:
            result = filenames[file_type]
        else:
            result = os.path.join(work_dir, jobid, filenames[file_type])
    return result

class Domain_hit:
    def __init__(self, domain_data):
        """
        Option <domain_data> could be either string or dictionary of previously obtained hit
        """
        if (type(domain_data) == type(str())) or (type(domain_data) == type(unicode())):
            #bZIP_Maf|179..269..4.3e-34..117.0
            self.domain_id = domain_data.split("|")[0]
            data = domain_data.split("|")[1].split("..")
            self.begin = int(data[0])
            self.end = int(data[1])
            self.evalue = float(data[2])
            self.score = float(data[3])
            if len(data) >= 6:
                self.hmm_begin = int(data[4])
                self.hmm_end = int(data[5])
            else:
                self.hmm_begin = None
                self.hmm_end = None
        if type(domain_data) == type(dict()):
            for key in domain_data.keys():
                setattr(self, key, domain_data[key])

    def get_plain_data_string(self):
        return "%i..%i..%s..%s" % (self.begin, self.end, self.evalue, self.score)

    def get_hit_length_in_protein(self):
        return self.end - self.begin + 1

class Protein:
    def __init__(self, domains_file_string = None, domain_as_dict = False):
        self.protein_id = "Unk"
        self.length = 0
        self.domain_hits = list()
        self.domain_formula = ""

        if domains_file_string != None:
            fields = domains_file_string.split("\t")
            #XP_008569256.1|Galeopterus_variegatus	291	bZIP_Maf|179..269..4.3e-34..117.0 Maf_N|86..119..1.2e-20..73.0	Maf_N bZIP_Maf
            self.protein_id = fields[0]
            self.length = int(fields[1])
            if len(fields) > 2: # Domains were found
                for domain_data in fields[2].split(" "):
                    if domain_as_dict:
                        self.domain_hits.append(Domain_hit(domain_data).__dict__)
                    else:
                        self.domain_hits.append(Domain_hit(domain_data))
                if domain_as_dict:
                    self.domain_hits = sorted(self.domain_hits, key = lambda k: k["begin"], reverse = False)
                else:
                    self.domain_hits = sorted(self.domain_hits, key = lambda k: k.begin, reverse = False)
                self.domain_formula = fields[3]

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def read_from_dict(self, dictionary):
        for key in dictionary.keys():
            if key == "domain_hits":
                domain_hits = list()
                for hit in dictionary[key]:
                    domain_hits.append(Domain_hit(hit))
                self.domain_hits = domain_hits
            else:
                setattr(self, key, dictionary[key])

def get_protein_data(domain_filename, as_dict = False):
    """
    Method returns (1) a list of <Protein> objects from the <domain_filename> (*.legend +
    *.domains output of the <udav_COGassign.py> script) and (2) strings for legend.
    If file was not found, returns None.
    If <as_dict> is 'True', returns a list of dicts rather than <Protein> objects.
    Can also take the file object as its first argument.
    """
    domain_file = None
    try:
        if not os.path.isfile(domain_filename):
            return None
        domain_file = open(domain_filename)
    except TypeError:  # This might be just file in memory
        domain_file = domain_filename
    legend_strings = list()
    reading_legend = False
    proteins = list()
    reading_proteins = False
    for string in domain_file:
        try:
            string = string.decode("utf-8")
        except:
            pass
        string = string.strip()
        if len(string) == 0:
            continue
        if re.match("#Legend", string): # Legend will be red
            reading_legend = True
            reading_proteins = False
        if re.match("#Output", string): # Main output will be red
            reading_proteins = True
            reading_legend = False
        if string[0] == "#":
            continue
        if reading_proteins:
            if as_dict:
                proteins.append(Protein(string, domain_as_dict = True).__dict__)
            else:
                proteins.append(Protein(string, domain_as_dict = False))
        if reading_legend:
            legend_strings.append(string.split("\t"))
    domain_file.close()
    return [proteins, legend_strings]

def get_hmmscan_command(work_dir, jobid):
    """
    Returns content of the file with the hmmscan string and None if no file was found
    """
    hmmscan_command_filename = get_filename(work_dir, jobid, "hmmscan_command")
    strings = list()
    try:
        hmmscan_command_file = open(hmmscan_command_filename)
        strings = hmmscan_command_file.readlines()
        hmmscan_command_file.close()
    except IOError:
        strings = None
    return strings

def get_all_proteins_from_request(work_dir, jobid):
    """
    Method returns a text version of a list of 'proteins' in FASTA-format corresponding to the 
    <work_dir> and <jobid>. Note that non-unique protein IDs were filtered out before.
    If no file was found, returns None.
    """
    fasta_filename = get_filename(work_dir, jobid, "fasta")
    result = ""
    try:
        strings = get_strings_from_fasta_file(fasta_filename)
        (proteins, non_unique) = udav_COGassign.get_proteins(strings, unique_required = False, quiet = True)
        fasta_strings = list()
        for p in proteins:
            fasta_strings.append(">%s\n" % p[0])
            fasta_strings.append("%s\n\n" % p[1])
        result = "".join(fasta_strings)
    except IOError:
        result = None
    return result

def get_strings_from_fasta_file(fasta_filename):
    """
    Returns: <strings> containing FASTA-formatted data from the <fasta_filename>
    """
    fasta_file = open(fasta_filename)
    curr_entry = "".join(fasta_file.readlines())
    fasta_file.close()
    if curr_entry.count(">") == 0:
        raise IOError("[..FATAL ERROR..] Input does not contain fasta-formatted sequences!")
    strings = curr_entry.split("\n")
    return strings

def check_and_write_fasta_file(fasta_data_string, work_dir, jobid):
    """
    This method takes <fasta_data_string> (input field of the user) and creates a file out of it
    in required <work_dir> with the <jobid>. Returns a name of this file and a dictionary of
    non-unique protein ids found in the request.
    """
    strings = fasta_data_string.split("\n")
    (proteins, non_unique) = udav_COGassign.get_proteins(strings, unique_required = True, quiet = True)
    non_unique_file = open(get_filename(work_dir, jobid, "non_unique"), "w")
    for n_u in non_unique.keys():
        non_unique_file.write("%s\n" % n_u)
    non_unique_file.close()

    fasta_filename = get_filename(work_dir, jobid, "fasta")
    fasta_file = open(fasta_filename, "w")
    id_file = open(get_filename(work_dir, jobid, "id_order"), "w")
    for p in proteins:
        id_file.write("%s\n" % p[0])
        fasta_file.write(">%s\n" % p[0])
        fasta_file.write("%s\n\n" % p[1])
    id_file.close()
    fasta_file.close()
    return (fasta_filename, non_unique)

def run_hmmscan_for_task(fasta_data_string, work_dir, jobid, hmmscan_path, profile_database_path):
    """
    This method will run an <hmmscan> for a given parameters in a background.
    <fasta_data_string>     - input field of the user (FASTA-formatted sequences in the 'Basic' format)
    <work_dir>              - main working directory
    <jobid>                 - current job id
    <hmmscan_path>          - path to the <hmmscan> binary file.
    <profile_database_path> - path to the profile database (check that it was pressed with the <hmmpress>)
    Returns: 1) the name of the FASTA-formatted file, 2) a dictionary of non-unique protein ids found in 
    the request. 3) domtblout output and 4) subprocess.Popen object
    """
    try:
        (temp_filename, non_unique) = check_and_write_fasta_file(fasta_data_string, work_dir, jobid)
    except:
        non_unique_file = open(get_filename(work_dir, jobid, "non_unique"), "w")
        non_unique_file.close()
        raise RuntimeError("Cannot write fasta file because of the broken sequence input!")
    domtblout_filename = get_filename(work_dir, jobid, "domtblout")
    output_filename = get_filename(work_dir, jobid, "out")
    hmmscan_stdout_filename = get_filename(work_dir, jobid, "hmmscan_stdout")
    hmmscan_stderr_filename = get_filename(work_dir, jobid, "hmmscan_stderr")
    args = [hmmscan_path, "--domtblout", domtblout_filename, "-o", output_filename, "--noali", profile_database_path, temp_filename]
    hmmscan_command_filename = get_filename(work_dir, jobid, "hmmscan_command")
    hmmscan_command_file = open(hmmscan_command_filename, "w")
    hmmscan_command_file.write(" ".join(args))
    hmmscan_command_file.close()
    with open(hmmscan_stdout_filename, "wb") as hmmscan_stdout, open(hmmscan_stderr_filename, "wb") as hmmscan_stderr:
        hmmscan_process = subprocess.Popen(args, stderr = hmmscan_stderr, stdout = hmmscan_stdout)
    return (temp_filename, non_unique, domtblout_filename, hmmscan_process)

def write_parameters_to_file(work_dir, jobid, params_domalyser):
    """
    Method writes all parameters in dictionary <params_domalyser> to a
    tab-separated file with required name for current <jobid>.
    """
    params_domalyser_file = open(get_filename(work_dir, jobid, "params_domalyser"), "w")
    params_domalyser_file.write("#Name\tValue\n")
    for key in params_domalyser.keys():
        params_domalyser_file.write("%s\t%s\n" % (key, params_domalyser[key]))
    params_domalyser_file.close()

def run_udav_tools(work_dir, jobid, database_filename, params_domalyser, fasta_data_string, hmmscan_path, profile_databases_path):
    """
    This method will create a working directory for a given <jobid>, write given
    <database_filename> and <params_domalyser> to specific files in this directory and
    run method <run_hmmscan_for_task>
    """
    curr_job_dir = os.path.join(work_dir, jobid)
    if os.path.isdir(curr_job_dir):
        raise IOError("Directory for the job '%s' already exists!" % curr_job_dir)
    os.mkdir(curr_job_dir)

    write_parameters_to_file(work_dir, jobid, params_domalyser)

    database_file = open(get_filename(work_dir, jobid, "database_name"), "w")
    database_file.write(database_filename)
    database_file.close()

    run_hmmscan_for_task(fasta_data_string, work_dir, jobid, hmmscan_path, os.path.join(profile_databases_path, database_filename))

def check_if_domtblout_file_is_ready(work_dir, jobid):
    """
    This method checks if the resulting domtblout output file of hmmscan in <work_dir> with <jobid> is ready.
    Returns a pair of <result> and <completion> percent. <result> is True if file is ready and properly finished,
    or False if not yet, but it is worth waiting (if no file was detected or file is empty, but errors are also empty)
    Raises error if there is no point in waiting any longer:
    <RuntimeError> if hmmscan gave an error into its STDERR
    """
    curr_job_dir = os.path.join(work_dir, jobid)
    if not os.path.isdir(curr_job_dir):
        raise IOError("Directory for the job '%s' is missing!" % curr_job_dir)
    domtblout_filename = get_filename(work_dir, jobid, "domtblout")
    result = None
    completion = 0 # estimated percent of task completion
    try:
        file_to_check = open(domtblout_filename)
        strings = file_to_check.readlines()
        file_to_check.close()
        result = False
        if strings[-1].strip() == "# [ok]":
            result = True
            completion = 100
        else:
            try:
                #COG1155              COG1155      588 Weevi_1212           -            528   4.5e-12   33.3   2.7   1   3       1.4       1.8   -5.0   5.4   514   573    28    87     5   151 0.77 Archaeal/vacuolar-type H+-ATPase catalytic subunit A/Vma1
                fields = re.split("\s+", strings[-2])
                query_id = fields[3]
                id_filename = get_filename(work_dir, jobid, "id_order")
                id_file = open(id_filename)
                strings = id_file.readlines()
                id_file.close()
                for i in range(len(strings)):
                    strings[i] = strings[i].strip()
                    if len(strings[i]) == 0:
                        continue
                    if strings[i] == query_id:
                        completion = int(100* i / len(strings))
                        break
            except IndexError:
                completion = None
    except IOError: # No file case
        result = False
        #print("I found no domtblout file yet!")
    except IndexError: # Empty file case
        result = False
        #print("Domtblout file appears to be empty")
    try:
        hmmscan_stderr_filename = get_filename(work_dir, jobid, "hmmscan_stderr")
        error_log_file = open(hmmscan_stderr_filename)
        error_strings = error_log_file.readlines()
        error_log_file.close()
        for string in error_strings:
            if re.match("Parse failed", string) or re.match("Error: File format problem, trying to open HMM", string) or re.match("Error: File existence", string) or re.match("Error: Sequence file", string):
                all_errors = "".join(error_strings)
                raise RuntimeError("%s" % all_errors)
    except IOError: # No hmmscan error log file case
        raise RuntimeError("No hmmscan error log file detected!")
    return (result, completion)

def get_database_name(work_dir, jobid):
    database_name_file = open(get_filename(work_dir, jobid, "database_name"))
    strings = database_name_file.readlines()
    database_name_file.close()
    return strings[0].strip()

def check_if_domains_file_is_ready(work_dir, jobid):
    """
    Checks if file <jobid>.domains exists in the corresponding job
    directory. Returns:
    None - if no file was detected
    True - if file exists and its last string corresponds to proper filling
    False - if file is empty or not yet filled until proper last string
    """
    result = None
    domains_filename = get_filename(work_dir, jobid, "domains")
    try:
        file_to_check = open(domains_filename)
        strings = file_to_check.readlines()
        file_to_check.close()
        if strings[-1].strip() == "# [ok]":
            result = True
        else:
            result = False
    except IOError: # No file case
        result = None
    except IndexError: # Empty file case
        result = False
    return result

def read_parameters_from_file(work_dir, jobid):
    """
    Method reads all parameters from a tab-separated file with required name
    for current <jobid> and returns them as a dictionary
    """
    params = dict()
    params_domalyser_file = open(get_filename(work_dir, jobid, "params_domalyser"))
    for string in params_domalyser_file:
        string = string.strip()
        if len(string) == 0:
            continue
        if string[0] == "#":
            continue
        fields = string.split("\t")
        params[fields[0]] = fields[1]
    params_domalyser_file.close()
    return params

def get_udav_state(work_dir, jobid):
    """
    This method checks if .domtblout and .domains files were obtained for a <jobid>
    and returns a <curr_state> dictionary with the following keys:
    "hmm" -> {
              "database_name",
              "state" ("finished"|"processing"|"error"),
              "errors" (None|string with hmmscan error)
              "output" (None|path to .domtblout file)
             }
    "da"  -> {
              "params" (dictionary of parameters),
              "state" ("finished"|"processing"|"error"),
              "errors" (None|dictionary of "parameter_name" -> "text"),
              "output" (None|dictionary with "domains_filename" -> path to .domains file,
                                             "domains_data" -> list returned by the <get_protein_data> method)
             }
    """
    curr_state = dict()
    curr_state["hmm"] = dict()
    try:
        curr_state["hmm"]["database_name"] = get_database_name(work_dir, jobid)
    except IOError:
        raise IOError("Job status is being checked before a job '%s' was started" % jobid)

    hmmscan_file_ready = None
    try:
        (hmmscan_file_ready, completion) = check_if_domtblout_file_is_ready(work_dir, jobid)
        if hmmscan_file_ready:
            curr_state["hmm"]["state"] = "finished"
            curr_state["hmm"]["errors"] = None
            curr_state["hmm"]["output"] = get_filename(work_dir, jobid, "domtblout", filename_only = True)
            curr_state["hmm"]["completion"] = completion
        else:
            curr_state["hmm"]["state"] = "processing"
            curr_state["hmm"]["errors"] = None
            curr_state["hmm"]["output"] = None
            curr_state["hmm"]["completion"] = completion
    except RuntimeError as e:
        curr_state["hmm"]["state"] = "error"
        curr_state["hmm"]["errors"] = {"hmm_params" : ["Программа hmmscan выдала следующие ошибки: '%s'" % str(e), "The hmmscan program encountered the following errors: '%s'" % str(e)]}
        curr_state["hmm"]["output"] = None
        curr_state["hmm"]["completion"] = 0

    seq_warnings = ""
    non_unique_file = open(get_filename(work_dir, jobid, "non_unique"))
    for string in non_unique_file:
        string = string.strip()
        if len(string) == 0:
            continue
        seq_warnings += "%s, " % string
    non_unique_file.close()
    if len(seq_warnings) != 0:
        seq_warnings = seq_warnings.strip(", ")
        curr_state["hmm"]["warnings"] = {"sequences" : ["Были обнаружены следующие не уникальные идентификаторы белков: %s. Все последовательности с одинаковыми идентификаторами (кроме первой) были удалены" % seq_warnings, "The following non-unique protein identifiers were detected: %s. All sequences with the same protein IDs except for the first occuring were removed" % seq_warnings]}
    else:
        curr_state["hmm"]["warnings"] = None

    curr_state["da"] = dict()
    p = read_parameters_from_file(work_dir, jobid)
    (p_c, errors, warnings) = check_web_params(p, as_str = True) # Checking, writing error messages
    curr_state["da"]["params"] = p_c
    curr_state["da"]["warnings"] = warnings
    if len(errors.keys()) != 0: # If any errors occur error messages are returned
        curr_state["da"]["state"] = "error"
        curr_state["da"]["errors"] = errors
        curr_state["da"]["output"] = None
    else:
        if hmmscan_file_ready == True:
            domains_ready = check_if_domains_file_is_ready(work_dir, jobid)
            if domains_ready == True:
                curr_state["da"]["state"] = "finished"
                curr_state["da"]["errors"] = None
                curr_state["da"]["output"] = dict()
                curr_state["da"]["output"]["domains_filename"] = get_filename(work_dir, jobid, "domains", filename_only = True)
                domains_filename = get_filename(work_dir, jobid, "domains")
                curr_state["da"]["output"]["domains_data"] = get_protein_data(domains_filename, True)
            else:
                if domains_ready == None: # No .domains file created yet
                    #write_parameters_to_file(work_dir, jobid, p_c)
                    run_domalyser_process(work_dir, jobid, p_c)
                curr_state["da"]["state"] = "processing"
                curr_state["da"]["errors"] = None
                curr_state["da"]["output"] = None
        else:
            if hmmscan_file_ready == False: # Worth waiting
                curr_state["da"]["state"] = "processing"
                curr_state["da"]["errors"] = None
                curr_state["da"]["output"] = None
            if hmmscan_file_ready == None: # HMMscan encountered fatal error
                curr_state["da"]["state"] = "error"
                curr_state["da"]["errors"] = None
                curr_state["da"]["output"] = None
    return curr_state

def run_domalyser_process(work_dir, jobid, p_c):
    """
    Method runs the <udav_COGassign_run.py> script with already checked <p_c>
    (dictionary of parameters)
    """
    script = os.path.join(p_c["udav_COGassign_run_path"], "udav_COGassign_run.py")
    args = [sys.executable, script, "--work_dir", work_dir, "--jobid", jobid, "--overlap", str(p_c["overlap"]), 
            "--evalue", str(p_c["evalue"]), "--evalue_type", str(p_c["evalue_type"])]
    if p_c["unite"] == "True":
        args.append("--unite")
    if "max_distance" in p_c:
        args.append("--max_distance")
        args.append(str(p_c["max_distance"]))
    if "max_hmm_overlap" in p_c:
        args.append("--max_hmm_overlap")
        args.append(str(p_c["max_hmm_overlap"]))
    domains_filename = get_filename(work_dir, jobid, "domains")
    if os.path.isfile(domains_filename): # Previous result is stored in the folder
        domains_file = open(domains_filename, "w")
        domains_file.close()
    domalyser_command_file = open(get_filename(work_dir, jobid, "domalyser_command"), "w")
    domalyser_command_file.write(" ".join(args))
    domalyser_command_file.close()
    domalyser_process = subprocess.Popen(args)

def check_and_run_domalyser_script(work_dir, jobid, params_domalyser):
    """
    Method checks if <udav_COGassign_run.py> script could be run, runs it if possible and
    returns a status of this run and a dictionary of errors. If the script started, parameters
    will be written to a corresponding file.
    (None, errors) - if parameters for run <params_domalyser> are wrong
    (None, dict()) - if hmmscan did not started working yet, so parameters check is irrelevant
    (False, dict()) - if .domains file is not yet ready but is in process
    (True, dict()) - if script started working
    """
    success = None
    errors = dict()
    if not os.path.isdir(os.path.join(work_dir, jobid)):
        raise IOError("No directory for jobid '%s' detected, cannot run domain_analyser" % jobid)
    else:
        (hmmscan_file_ready, completion) = check_if_domtblout_file_is_ready(work_dir, jobid)
        if hmmscan_file_ready:
            domains_ready = check_if_domains_file_is_ready(work_dir, jobid)
            if domains_ready == False: # This means that file was created but it is not yet filled, the process is running
                success = False
                raise RuntimeError("Cannot run domain_analyser because .domains file is not yet finished!")
            else:
                (p_c, errors, warnings) = check_web_params(params_domalyser, as_str = True) # Checking, writing error messages
                write_parameters_to_file(work_dir, jobid, params_domalyser)
                if len(errors.keys()) != 0:
                    domains_filename = get_filename(work_dir, jobid, "domains")
                    if os.path.isfile(domains_filename):
                        os.remove(domains_filename)
                    success = None
                else:
                    run_domalyser_process(work_dir, jobid, p_c)
                    success = True
    return (success, errors)

def check_single_parameter(web_dict, key, default, method_to_transform, error_messages, errors, warnings, min_value = None, max_value = None, as_str = False):
    value = None
    if key in web_dict:
       try:
           value = method_to_transform(web_dict[key])
           if min_value != None:
               if value < min_value:
                   warnings[key] = ["Поданная величина %s меньше, чем минимально допустимая (%s)" % (value, min_value), "Provided value %s is less than minimal allowed (%s)" % (value, min_value)]
                   value = min_value
           if max_value != None:
               if value > max_value:
                   warnings[key] = ["Поданная величина %s больше, чем максимально допустимая (%s)" % (value, max_value), "Provided value %s is larger than maximum allowed (%s)" % (value, max_value)]
                   value = max_value
       except ValueError:
           value = web_dict[key]
           errors[key] = [error_messages[0] % web_dict[key], error_messages[1] % web_dict[key]]
    else:
        value = method_to_transform(default)
        warnings[key] = ["Значение не было подано, установлено значение по умолчанию %s" % default, "The value was not provided, default value %s is used" % default]

    if as_str:
        value = str(value)
    return value

def check_web_params(web_dict, as_str = False):
    """
    Method will check if parameters passed by the <web_dict> are correct.
    If <as_str> is True, parameters values will be returned as strings.
    Also returns dictionaries of <errors> and <warnings>.
    """
    new_dict = dict()
    errors = dict()
    warnings = dict()
    new_dict["evalue"] = check_single_parameter(web_dict, "evalue", "1e-05", float, ["Поданное значение '%s' не является дробным числом или числом с плавающей точкой", "Non-float value '%s' was given"], errors, warnings, min_value = 1e-190, max_value = 10, as_str = as_str)
    new_dict["overlap"] = check_single_parameter(web_dict, "overlap", 20, int, ["Поданное значение '%s' не является целым числом", "Non-int value '%s' was given"], errors, warnings, min_value = 0, max_value = 100, as_str = as_str)

    if "evalue_type" in web_dict:
        if not web_dict["evalue_type"] in ["i", "c"]:
            new_dict["evalue_type"] = "i"
        else:
            new_dict["evalue_type"] = web_dict["evalue_type"]
    else:
        new_dict["evalue_type"] = "i"

    if "unite" in web_dict:
        if not str(web_dict["unite"]) in ["True", "False"]:
            new_dict["unite"] = "False"
        else:
            new_dict["unite"] = web_dict["unite"]
            if new_dict["unite"] == "True":
                new_dict["max_distance"] = check_single_parameter(web_dict, "max_distance", 50, int, ["Поданное значение '%s' не является целым числом", "Non-int value '%s' was given"], errors, warnings, min_value = 0, as_str = as_str)
                new_dict["max_hmm_overlap"] = check_single_parameter(web_dict, "max_hmm_overlap", 20, int, ["Поданное значение '%s' не является целым числом", "Non-int value '%s' was given"], errors, warnings, min_value = 0, max_value = 100, as_str = as_str)
    else:
        new_dict["unite"] = "False"

    try:
        script = os.path.join(web_dict["udav_COGassign_run_path"], "udav_COGassign_run.py")
        if not os.path.isfile(script):
            raise IOError("No script file found at '%s'!" % script)
        new_dict["udav_COGassign_run_path"] = web_dict["udav_COGassign_run_path"]
    except IOError:
        errors["udav_COGassign_run_path"] = ["По данному пути '%s' не найдено файла <udav_COGassign_run.py>" % web_dict["udav_COGassign_run_path"], "No file <udav_COGassign_run.py> found under path '%s'" % web_dict["udav_COGassign_run_path"]]
    except KeyError:
        errors["udav_COGassign_run_path"] = ["Не указан путь к файлу со скриптом <udav_COGassign_run.py>", "Path to file with the <udav_COGassign_run.py> script is not specified"]
    return (new_dict, errors, warnings)

def get_hmm_database_list(top_hmm_directory_path):
    """
    Gets a list of *.hmm inside the <top_hmm_directory_path> together with information
    from the '_HMM_databases.txt' file which should also be placed in the same directory.
    Returns: list of [path_to_hmm_database, [ru_hmm_database_descr, en_hmm_database_descr], 
    hmm_database_human_name],
    if not information is found, two latter will be None.
    """
    database_list_filename = os.path.join(top_hmm_directory_path, "_HMM_databases.txt")
    known_databases = dict()
    database_order = list()
    try:
        database_list_file = codecs.open(database_list_filename, encoding = "utf-8")
        for string in database_list_file:
            string = string.strip()
            if len(string) == 0:
                continue
            if string[0] == "#":
                continue
            fields = string.split("\t")
            curr_filename = fields.pop(0)
            known_databases[curr_filename] = [fields[0].split("|"), fields[1]]
            database_order.append(curr_filename)
        database_list_file.close()
    except IOError:
        raise IOError("FATAL ERROR: file '%s' not found!" % database_list_filename)
    dir_content = os.listdir(top_hmm_directory_path)
    #print (known_databases)
    found_databases = dict()
    for d in dir_content:
        if re.search("\.hmm$", d) != None: # This is an *.hmm file
            found_databases[d] = True

    result = list()
    for filename in database_order:
        if filename in found_databases:
            result.append([filename] + known_databases[filename])
            found_databases.pop(filename)
        else:
            print ("WARNING: database %s is stated in the file, but not found!" % filename)
    for filename in found_databases.keys():
        result.append([filename, None, None])
    return result

def get_curr_domains_file(path_to_domains_file, options):
    """
    Provide a path to .domains file and <options> containing a 'colormap' key: a dictionary 
    linking domain name and its current color (will replace legend strings in the .domains file)

    Returns: string with all data for the new .domains file
    """
    # 1) Reading <proteins> and <legend_strings>
    proteins = None
    legend_strings = None
    if os.path.isfile(path_to_domains_file):
        (proteins, legend_strings) = get_protein_data(path_to_domains_file)
    else:
        raise IOError("The following path is not a file: '%s'\n" % path_to_domains_file)
    # 2) Extracting actual <domain_to_color> data from the given <options['colormap']>
    domain_to_color = dict()
    try:
        for single_domain_data in options["colormap"].split("\n"):
            domain_data = single_domain_data.split(" ")[0]
            color = single_domain_data.split(" ")[1]
            domain_to_color[domain_data] = color
    except KeyError:
        raise KeyError("Bad options passed; they must contain key 'colormap': %s" % options)
    # 3) Changing <legend_strings>
    for i in range(len(legend_strings)):
        try:
           legend_strings[i][2] = domain_to_color[legend_strings[i][1]]
        except KeyError:
           raise KeyError("Domain '%s' is missing in the provided colormap!" % legend_strings[i][1])

    # 4) Getting result string
    result_strings = list()
    curr_domains_file = open(path_to_domains_file)
    reading_legend = False
    for string in curr_domains_file:
        if re.match("\#Number", string) != None: # First string of the legend
            reading_legend = True
            result_strings.append(string)
            for legend_string in legend_strings:
                result_strings.append("%s\n" % "\t".join(legend_string))
            continue
        if reading_legend and (string[0] == "#"): # Legend finished
            reading_legend = False
        if not reading_legend:
            result_strings.append(string)
         
    return "".join(result_strings)

def make_svg(source, options, as_string = False):
    """
    Writes a custom svg file from the specified data. If <source> is a list, it should
    contain a <work_dir> and <jobid> whire a proper .domains file should reside. Or provide
    a .domains filename

    Parameters are specified with <options> dict with the following keys:
    * 'colormap' - dictionary linking domain name and its current color (will replace legend
    strings in the .domains file or given data).
    * 'details' - should be either an empty string or 'evalue' or 'score'
    * 'align' - should be ether an empty string or name of the domain to align figure to
    * 'as_scheme' - if True, result will contain schematic view of each protein domains

    If <as_string> is True, no file will be created, but the result is returned as an svg string.
    """
    proteins = None
    legend_strings = None
    if type(source) == type(list()): # Input in the form of proteins and strings (no job was run on the server)
        try:
            proteins = list()
            for p in source[0]:
                new_prot = Protein()
                new_prot.read_from_dict(p)
                proteins.append(new_prot)
            legend_strings = source[1]
        except IndexError:
            raise IndexError("Wrong <source> given to the <make_svg()> method: '%s'" % source)
        if as_string == False:
            raise IOError("File creation required from the <make_svg()> method, but no directory to put it is provided!")
    elif type(source) == type(str()): # Input is supposed to be a filename
        if os.path.isfile(source):
            (proteins, legend_strings) = get_protein_data(source)
        else:
            raise IOError("The following path is not a file: '%s'\n" % source)
    add_measure = None
    domain = None
    natural = True
    domain_to_color = dict()
    try:
        for single_domain_data in options["colormap"].split("\n"):
            domain_data = single_domain_data.split(" ")[0]
            color = single_domain_data.split(" ")[1]
            domain_to_color[domain_data] = color
        if options["details"] != "":
            add_measure = options["details"]
        if options["align"] != "":
            domain = options["align"]
        if "as_scheme" in options:
            natural = not options["as_scheme"]
    except KeyError:
        raise KeyError("Bad options passed; they must contain keys 'colormap', 'details', 'align': %s" % options)

    for i in range(len(legend_strings)):
        try:
           legend_strings[i][2] = domain_to_color[legend_strings[i][1]]
        except KeyError:
           raise KeyError("Domain '%s' is missing in the provided colormap!" % legend_strings[i][1])
        legend_strings[i] = "\t".join(legend_strings[i])
    proteins_new = list()
    for p in proteins:
        pseudosequence = "A" * p.length
        pseudofasta_string = ""
        new_domain_data = list()
        p.domain_hits.sort(key = lambda k: k.get_hit_length_in_protein(), reverse = True)
        for hit in p.domain_hits:
            new_domain_data.append((hit.domain_id, hit.get_plain_data_string()))
        proteins_new.append([p.protein_id, pseudosequence, new_domain_data, pseudofasta_string])
    svg_string = None
    if as_string:
        svg_string = udav_COGassign.write_svg_scheme(None, proteins_new, domain_to_color, legend_strings, domain = domain, natural = natural, add_measure = add_measure)
    else:
        jobid = os.path.basename(os.path.dirname(source))
        work_dir = os.path.dirname(os.path.dirname(source))
        output_filename = get_filename(work_dir, jobid, "svg")
        svg_string = udav_COGassign.write_svg_scheme(output_filename, proteins_new, domain_to_color, legend_strings, domain = domain, natural = natural, add_measure = add_measure)
    return svg_string

if __name__ == '__main__':
    print ("Module <udav_COGassign_utils> (version %s) was executed as a script." % curr_version)
    print ("It will check if it has access to data and if your options are valid.")
    if len(sys.argv) > 2:
        print ("FATAL ERROR: wrong number of command line arguments!")
        print ("Please provide either no arguments or 1 command line argument:")
        print ("    1) Path to the *.domains output of the <udav_COGassign.py>.")
        sys.exit()
    domains_filename = None
    if len(sys.argv) == 2:
        domains_filename = sys.argv[1]
    else:
        #prot1 = ">Dole_0603\nMSDNIGKVVQVMGPVVDVEFEPGKLPAILTALLITNTVINDEADNLVVEVAQHLGDNVVRTIAMDVTDGLVRGMPVKDTGAPITMPVGAASLGRVLNVVGKPVDGLGPVSREKTMPIHRPAPLFTEQDTSVNVLETGIKVIDLLVPFPRGGKMGLFGGAGVGKTVIMMEMVNNIAMQHGGISVFAGVGERTREGNDLYHEMKDSGVLPKAALIYGQMTEPPGARARVALSALTCAEYFRDVEGQDVLIFIDNIFRFTQAGAEVSALLGRIPSAVGYQPTLAVDLGGLQERITSTDKGSITAVQCVYVPADDLTDPAPATTFAHLDGTVVLSRQIAELGIYPSVDPLDSTSRILDAAYIGEEHYRVAREVQQTLQKYKELQDIIAILGMDELSDEDKVTVERARKLQRFLSQPFHVAEVFTGKPGSYVKIEDTVRSFKEICDGKHDDLPESAFYMVGSIEEAVAKAKG\n"
        #prot2 = ">Theam_1659\nMQIRAEEISELIRKQIEEFEASVNLDETGIVIKVGDGVARVYGLENVEYGEVVEFEDGTEGVAFNLEEDNVGVVLLGEGRGIVEGGKAKRTGRILDMPVGDGLIGRVLDPLGNPIDGKGDIEYTERRAVERIAPGIVTRKPVHEPLQTGIKAIDALIPIGRGQRELIIGDRQTGKTTVAIDTILNQKREGVICVYCAVGQKRSTVAQTIQLLKELGAMDYTIVISATASDPAALQYLAPYAACTVAEYFRDTGRAALIVYDDLSKQAVAYREMSLLLRRPPGREAYPGDVFYLHSRLLERAAKLNDELGAGSLTALPIVETKAGDISAYIPTNVISITDGQIFLETDLFYKGQRPAINVGLSVSRVGGAAQIKAMKQVAGKLRLELARYRELEAFAQFASDLDPATRAQLERGRRLMELLKQPPHKPIPVEKQIVAFFAAINGYLDDIPVEAVTKFEWELYAFMDAKHPEILKEILEKKKLDDELTKKLHEAIKEFKATFTA\n"
        #prot3 = ">plu3767\nMKLSLDHIPGKMRHAINECRLIQIRGRVTQVTGTLLKAVIPGVRIGELCHLRNPDNTLSLLAEVIGFQQHQALLTPLGEMFGISSNTEVSPTGAMHQVGVGDYLLGQVLDGLGNPFSGGQLPEPQAWYPVYRDAPAPMSRKRIEHPLSLGVRAIDGLLTCGEGQRMGIFAAAGGGKSTLLSTLIRSAEVDVTVLALIGERGREVREFIESDLGEEGLKRSVLVVATSDRPAMERAKAGFVATSIAEYFRDQGKRVLLLMDSVTRFARAQREIGLAAGEPPTRRGYPPSVFAALPRLMERAGQSDKGSITALYTVLVEGDDMTEPVADETRSILDGHIILSRKLAAANHYPAIDVLRSASRVMNQIITPEHQAQAGLLRKWLAKYEEVELLLQIGEYQKGQDPVADNAIAHIEAIRNWLRQGTHEPSDLPQTLAQLQQITK\n"
        #protein_strings = prot1 + prot2 + prot3
        dummy_file = open("dummy.fasta")
        protein_strings = "".join(dummy_file.readlines())
        dummy_file.close()
        work_dir = "D:\\UdavBackup\\_Complete_genomes\\_Data\\COG2020_analysis\\_domain_analyser"
        jobid = "lalala"
        hmmscan_path = os.path.join(work_dir, "hmmscan.exe")
        params_domalyser = {"evalue" : "1e-5", "overlap" : "-1", "max_distance" : "50", "max_hmm_overlap" : "20", "evalue_type" : "i", "unite" : "True", "udav_COGassign_run_path" : "D:\\UdavBackup\\_Complete_genomes\\_scripts\\special"}
        run_udav_tools(work_dir, jobid, "_dummy.hmm", params_domalyser, protein_strings, hmmscan_path, work_dir)
        import time
        for i in range(5):
            params = get_udav_state(work_dir, jobid)
            #print ("%i half-seconds, hmmscan completion = %i, domalyser run: '%s'" % (i, completion, success))
            print ("%i half-seconds, hmmscan completion = %i" % (i, params["hmm"]["completion"]))
            #if len(errors.keys()) != 0:
            #    print (errors)
            time.sleep(0.5)
        params_domalyser = {"evalue" : "1e-5", "overlap" : "-1", "max_distance" : "50", "max_hmm_overlap" : "20", "evalue_type" : "i", "unite" : "True", "udav_COGassign_run_path" : "D:\\UdavBackup\\_Complete_genomes\\_scripts\\special"}
        #(success, errors) = check_and_run_domalyser_script(work_dir, jobid, params_domalyser)
        params = get_udav_state(work_dir, jobid)
        print (params["hmm"])
        print (params["da"]["state"])
        print (params["da"]["params"])
        print (params["da"]["errors"])
        print (params["da"]["warnings"])
        #print (params["da"]["output"]["domains_filename"])
        #params_domalyser = {"evalue" : "1e-10", "overlap" : "90", "max_distance" : "50", "max_hmm_overlap" : "20", "evalue_type" : "i", "unite" : "True", "udav_COGassign_run_path" : "D:\\UdavBackup\\_Complete_genomes\\_scripts\\special"}
        #status = check_and_run_domalyser_script(work_dir, jobid, params_domalyser)
        #print (status)
        #time.sleep(2.5)
        #params = get_udav_state(work_dir, jobid)
        #print (params)
        #options = {"colormap":"COG1155 #567834\nCOG1156 #111111\nCOG1157 #880088\nCOG0055 #008888\nCOG0056 #998877", "details":"score", "align":"","as_scheme":False}
        #result = make_svg(params["da"]["output"]["domains_filename"], options)
    hmm = get_hmm_database_list("D:\\UdavBackup\\_Complete_genomes\\_Data\\COG2020_analysis\\_domain_analyser")    
    print (hmm)
    options = dict()
    options["colormap"] = "COG1157 #71C837\nCOG0055 #888888\nCOG0056 #FF00FF\nCOG1155 #18B7C8"
    result = get_curr_domains_file(domains_filename, options)
    myoutput = open("file.txt", "w")
    myoutput.write(result)
    myoutput.close()        