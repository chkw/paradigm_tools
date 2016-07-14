# paradigm_tools
tools relating to preparing data for Paradigm and post-processing Paradigm results
Running PARADIGM via five3 Website
===

PARADIGM is a factor graph method for inferring biological pathway activity from gene expression and gene copy number sample data. More information about PARADIGM is available at:

  - <http://bioinformatics.oxfordjournals.org/content/26/12/i237>
  - <http://bioinformatics.oxfordjournals.org/content/29/13/i62>

data requirements
---

  - pathway file
    - Here at UCSC, we use our superpathway. The current version is `superpathway_4.0.no_chemicals.tab`. You could create your own pathway file by following the format described at <https://dna.five3genomics.com/paradigm/help>.
  - evidence files (expression data & copy number data)
    - The sample data are provided in tab-delimited 2-D matrix format.
	  - samples in the columns
	  - genes in the rows
	- The gene names must match the gene names in the pathway file (usually HUGO symbol). 
	- The sample names must be common across data files as well as in the same order.
	  - use `commonColumnOrder.py`
	- **Are the features (genes) required to be common across data files as well as in the same order?**
	  - use `commonRowOrder.py`
	- normalizing expression data
	  - If no match normal data is available, normalize against the cohort itself using a gene-wise, mean-center method.
	    - `Rscript normalize.r mRNA_matrix.tab mean_shift row > mRNA_matrix.meanshifted.tab`
	  - If normal data is available, then normalize against the normal data.
	    - **(Do we have a script for doing this?)**
    - **Is there any normalization requirement for copy number data?**

parameters
---

  - **starting parameters and training parameters?**

submitting the job
---

Submit the job via website at <http://paradigm.five3genomics.com/>.

PARADIGM results
---

(What to do with the results?)

(modulation filter)

(constituent pathways most enriched for changed activity level)

