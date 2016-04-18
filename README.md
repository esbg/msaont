#MSAOnt Population Software
------
This software generates an instance of the multiple sequence alignment ontology
(MSAOnt) from an aligned fasta or cma file.

Requires Python version 2.7

License: MIT

#Python dependencies
------
Please ensure the following Python packages are installed via ``pip``:
- `Biopython <http://biopython.org/wiki/Main_Page>`
- `RDFLib <https://github.com/RDFLib/rdflib>`
- `progress <https://pypi.python.org/pypi/progress>`
- `biocma <https://github.com/etal/biocma>`

#Note
------
If you are populating MSAOnt with a cma file, the first sequence must be an 
aligned consensus (with no deletions).

