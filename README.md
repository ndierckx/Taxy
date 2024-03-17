# Taxy
Taxonomy assignment of ASVs

## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/Taxy/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

<a href="https://groups.oist.jp/macc" target="_blank"><img border="0" src="https://pbs.twimg.com/profile_images/1172654211579990018/LuQCVXkn_400x400.jpg" width=auto height="95" ></a> 
<a href="https://groups.oist.jp/macc" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/e/ef/OIST_logo.png" width=auto height="95" ></a> 

## Instructions

### 1. Install dependencies

- Install BLAST
- Install MAFFT
- Install Perl modules: MCE::Child && MCE::Channel

#### With Conda:

<code>conda create -n taxy -c conda-forge -c bioconda perl</code>

<code>conda install blast</code>

<code>conda install mafft</code>

<code>conda install perl-mce</code>
  
## Configuration file

This is an example of a configuration file for the assembly of a chloroplast.
To make the assembler work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

**1. Example of configuration file:**
<pre>
Project name              = Test
Combined reads or ASVs    = /path/to/reads/or/ASVs/ASVs.fasta
Forward reads             = /path/to/reads/reads_1.fastq (at the moment only the ASV option is available)
Reverse reads             = /path/to/reads/reads_2.fastq (at the moment only the ASV option is available)
Keep read ids             = 
Nucleotide database       = /path/to/nucleotide/database/from/NCBI/nt
Taxonomy database         = /path/to/taxanomy/database/from/NCBI/taxonomy.xml
Taxonomy only             = yes
Threads                   = 4
Output path               = /path/to/output/folder/
</pre>

**2. Explanation parameters:**
<pre>
#Project name              = Choose a name for your project, it will be used for the output files.
#Combined reads or ASVs    = /home/nicolas/Perl/OIST/eDNA/Michael/NC_asv_table.fasta
#Forward reads             = The path to the file that contains the forward reads (not necessary when there is a Combined or ASV file)
#Reverse reads             = The path to the file that contains the reverse reads (not necessary when there is a Combined or ASV file)
#Keep read ids             = When yes, the read ids from the fasta or fastq files are used in the output files, otherwise they will be changed to numbers (yes/no)
#Nucleotide database       = /path/to/nucleotide/database/from/NCBI/nt (https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl) (Other databeses with the same structure can also be used)
#Taxonomy database         = /home/nicolas/Perl/OIST/eDNA/taxonomy.xml (http://ftp.ebi.ac.uk/pub/databases/ena/taxonomy/taxonomy.xml.gz)
#Taxonomy only             = When yes, ASVs will directly be used for taxonomy assignment without prior clustering. (yes/no)
#Threads                   = Increasing the number of cores will speed up the runtime.
#Output path               = You can change the directory where all the output files wil be stored.
</pre>
</html>
