This repository contains pymol scripts received from Robert Campbell of Queen's University on 1/31/2024 for public usage.  Note that each script contains a copyright declaration.  

These scripts are those described on a now lost Queen's University webpage.  A shallow copy of the webpage exists on the wayback machine: https://web.archive.org/web/20160104232715/http:/pldserver1.biochem.queensu.ca/~rlc/work/pymol

A copy of the table from the website is available in this repo and called 'robert_campbell_script_summary.xlsx'

Robert says "I haven't tested all of these recently so I would not want to verify that they all work as originally designed, but most should.  The ones I'm less sure of are the ones in the cctbx subdirectory as those clearly have not been updated to python3 and I haven't had cctbx installed for a long time!"

from the original website:

All of these scripts require loading into PyMOL before use. You do this with the run command:
```
run scriptname.py
```
Then you can use the command according to the instructions provided (see 'robert_campbell_script_summary.xlsx' file).
Some of these require cctbx the Computational Crystallography Toolbox (marked with cctbx) and NumPy (Numerical Python) (numpy) and some require PyMOL version 0.8 (or higher) due to the use of the pymol.vfont module (PyMOL +v0.80). The following have only really been tested on the SVN versions of PyMOL.
