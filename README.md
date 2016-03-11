Design Guide RNA - CSC8311 Assignment (20015/2016) README file
==============================================================

This tool for the CSC8311 Assignment is developed by Django with Python and 
Biopython. 

It tries to find guide RNA in a sequence similar to the target sequence in a 
genome and stores the results in database. 

Firstly, the query sequence (genome of a species where we want to find gRNAs) 
must be selected from the already stored sequences (currently one). 
After that, the target sequence (gene or DNA sequence) must be applied to find 
homologous sequence in the genome using BLAT standalone command. 
If homologous sequence(s) is/are found (hit(s)) in the genome, then the tool 
finds the prospective guide RNAs in the hits based on the PAM sequence 
selected on the main page. 

Currenlty, only the wild-type *blunt* Cas9 is supported.

Requirements
============

The Design gRNA Tool (grna) is currently supported and tested on Ubuntu LTS 
14.04 (x64) Linux and Mac OS X > 10.10. 
The other requirements are the following:
 
- Django 1.9.3 -- see https://www.djangoproject.com

  The primary develompent platform for Design gRNA tool.
  
- Python 2.7.x -- see http://www.python.org

  Django requires python for development.
  
- Biopython 1.65 -- see http://biopython.org/wiki/Main_Page

  Biological tool written in Python and used for the assignment.

- Python new regex module - see https://pypi.python.org/pypi/regex

  To install run _sudo pip install regex_



Dependencies
============

- blat - Standalone BLAT v. 36x1, see https://genome.ucsc.edu/FAQ/FAQblat.html

  Blat is fast sequence search command line tool for Linux x64 and Mac OS X 
  and is bundled under the ./utils directory and it is free for academic or 
  education purposes.
  
- git - standalone for installing Design gRNA by cloning it from github, see 
https://git-scm.com/downloads

- pip - Django and BioPython requires _pip_ to be installed.
 
- python-dev - Biopython requires _python-dev_ package to be installed on Ubuntu

- regex - gRNA tool requires the new regular expression library _regex_.

Installation and use
======================

Firts, **make sure that proper version of Python, Django, git standalone and 
Biopython are 
installed correctly**. Then, installation can be done by unpakcaging  
gzipped tar or clone from github. See example below to install grna 
and its requirements in Ubuntu.

Install requirements and dependencies:

        $ python --version
        Python 2.7.6
        
        $ sudo apt-get update
        
        $ sudo apt-get install git
              
        $ sudo apt-get install python-pip
       
        $ sudo apt-get install python-dev 
        
        $ sudo pip install Django==1.9.3
        Downloading/unpacking Django==1.9.3
       
        $ sudo pip install Biopython==1.65
        Downloading/unpacking Biopython==1.65
        
        $ sudo pip install regex

Installation from github:
  
       git clone https://github.com/ilap/CSC8311A1


Installation from the package:
       
        tar zxvf CSC8311A1.tgz
        cd CSC8311A1

Run the site:

       cd CSC8311A1
       python manage.py runserver
            Starting development server at **http://127.0.0.1:8000/**
       
Access it using browser opening the IP:port showed above, and then proceed 
the following steps below:

1. Select the species genome (currently only _B. Subtillis strain 168_),

2. Insert the target sequence to the Species text area to find homologous sequence in the selected genome.
   For testing, copy and paste the contain of the  *amyE* gene of *B. 
subtillis* fasta file (amyE_B_Subtilis_strain_168.fa).
3. Choose the up/down stream offset for extending the search range by the offset in the genome.
   The 0 means, that the search will be run only on the exact length of the found target sequence.
   The deafult (2000) means that the search will start 2000 nucleotides upstream and finish at 2000 nucleotides
   downstream of the selected gene in the genome.
4. Select PAM sequence. Currently, just some of the wild-type Cas9's PAMs (NGG and NAG) are added to the database.
5. And finally click the **__Search gRNA__** button.

The result (see screenshots) should be the list of the found guide RNAs in the 
genome, based on the target sequence. There is **NO** any statistical analysis in 
this base assignment.

See, screenshot below:

![Screenshot](./misc/Screenshot.png?raw=true "Screenshot of Design gRNA Tool")
![Screenshot1](./misc/Screenshot1.png?raw=true "Result's screenshot of Design gRNA Tool")

Testing
=======

Run the following command to run test unit for grna.

      python manage.py  test
      Creating test database for alias 'default'...
      ...
      ----------------------------------------------------------------------
      Ran 3 tests in 0.006s

      OK
      Destroying test database for alias 'default'...

Distribution Structure
======================

- README       -- This file.
- mysite/      -- The Django ROOT site.
- grna/        -- The Design gRNA (grna) tool site.
- db.sqlite3   -- The databese for the grna.
- manage.py    -- The control file of Django.
- sequences/   -- The miscellaneous genome and target sequences used for 
  developing and testing grna.
- utils/       -- The required BLAT standalone program for Linux x64 and Mac 
OS X.
- misc/        -- Miscellaneous stuffs e.g. screenshots.
