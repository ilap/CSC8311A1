Design Guide RNA - CSC8311 Assignment (20015/2016) README file
==============================================================


This tool for the CSC8311 Assignment is developed by Django with Python and 
Biopython. 

It try to find guide RNA in a sequence similar to the target sequence in a 
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
  
- Python 2.7 -- see http://www.python.org

  Django requires python for development.
  
- Biopython 1.65 -- see http://biopython.org/wiki/Main_Page

  Biological tool written in Python and used for the assignment.



Dependencies
============

- blat - Standalone BLAT v. 36x1, see https://genome.ucsc.edu/FAQ/FAQblat.html

  Blat is fast sequence search command line tool for Linux x64 and Mac OS X 
  and is bundled under the ./utils directory and it is free for academic or 
  education purposes.
  
- git - standalone for installing Design gRNA by cloning it from github, see 
https://git-scm.com/downloads

Installation and use
======================

Firts, **make sure that proper version of Python, Django, git standalone and 
Biopython are 
installed correctly**. Then, installation can be done by unpakcaging  
gzipped tar or clone from github.

Installation from github:
  
       git clone https://github.com/ilap/CSC8311A1


Installation from the package:
       
        tar zxvf CSC8311A1.tgz
        cd CSC8311A1

Run the site:

       cd CSC8311A1
       python manage.py runserver
            Starting development server at **http://127.0.0.1:8000/**
       
Access it using browser opening the IP:port showed above.

-- Select the species genome,
-- Insert the target sequence to find honologous sequence in the selected 
genome.
-- Choose the up/down stream offset for extending the search range by 
the offset in the genome. The 0 means, that the search will be run only on the 
exact length of the found target sequence.
-- Seelct PAM
-- And finally clikc the *Search gRNA* button.

See, screenshot below:

![Screenshot](./misc/Screenshot.png?raw=true "Screnshot of Design gRNA Tool")

Testing
=======

**Not yet implemented**


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
