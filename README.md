# COGassign
This is a basic script behind the web-version of the DomainAnalyser tool (https://boabio.belozersky.msu.ru/DomainAnalyser) which can be used
for local searches with any custom database and a set a sequences.

# Data preparatopn
* Please provide sequences in fasta format. **Sequence IDs should be unique** in a single run, or HMMer will fail to run them.
* Use ```hmmpress``` on the profile HMM database. Profile HMM database can be obtained, e.g., from Pfam (https://www.ebi.ac.uk/interpro/download/pfam/) or for COGs (https://boabio.belozersky.msu.ru/tools).
* You should install ```hmmer``` for your system (http://hmmer.org/). Windows release is available here: http://eddylab.org/software/hmmer3/3.1b2/ 
(if it fails to run, consider putting ```cygwin1.dll``` into a directory with binary files).
* (OPTIONAL) Install ```svgwrite``` Python module, or no SVG output will be created.

# Common usage (Windows):
This run would use ```mydata.fasta``` sequence file, search for Pfam domain hits in it and: 
* remove hits which are worse than ```0.01``` by e-value;
* filter them on ```5%``` overlap threshold (see https://boabio.belozersky.msu.ru/DomainAnalyser for options description);
* try to unite domain parts (under default options of 'Max Distance' and 'Max HMM Overlap (%)');
* domain scheme (with domain sizes proportional to portions of a protein occupied by their hit) will be written to ```mydata.domains.svg``` and sorted with the first occurence of domain ```DNA_pol_A```;
* use color file ```mydata.colors``` to color domains in the resulting SVG (color file is a simple tab-delimited two-column text file where first columns is domain name, e.g. ```DNA_pol_A```, and second column is HEX-coded color).
```
set hmmscan_path=D:\bioinf_soft\hmmer-3.1b2-cygwin64\binaries
set hmm_database=D:\Databases\HMM_database\Pfam-A_37.2.hmm
set prefix=mydata
udav_COGassign.py -p %hmm_database% -i %prefix%.fasta -o %prefix%.domains -l %prefix%.legend -e 0.01 -f 5 -c %prefix%.colors --unite --hmmscan %hmmscan_path% -d "DNA_pol_A"
```
Run the script with the ```-h``` option to get a complete list of available options.

! Please note that the script pretends to be clever and checks if ```%prefix%.table.domains``` exists. If the file is there, HMMer will not run, instead this file will be used for further work.
