# ENCODE Imputation Challenge Scoring and Validation Scripts

## Installation

1) [Install Conda](https://docs.conda.io/en/latest/miniconda.html) first.

2) Install `numpy`, `scikit-learn` and `pyBigWig`.
	```bash
	$ conda install -c bioconda numpy scikit-learn pyBigWig sqlite scipy
	```

## Example (hg38)

1) Download ENCFF622DXZ and ENCFF074VQD from ENCODE portal.
	```bash
	$ mkdir -p test/hg38 && cd test/hg38
	$ wget https://www.encodeproject.org/files/ENCFF622DXZ/@@download/ENCFF622DXZ.bigWig
	$ wget https://www.encodeproject.org/files/ENCFF074VQD/@@download/ENCFF074VQD.bigWig
	```

2) Convert it to numpy array.
	```bash
	$ python bw_to_npy.py test/hg38/ENCFF622DXZ.bigWig
	$ python bw_to_npy.py test/hg38/ENCFF074VQD.bigWig
	```

3) Run it. If you score without a variance `.npy` file specified as `--var-npy`, then `msevar` metric will be `0.0`.
	```bash
	$ python score.py test/hg38/ENCFF622DXZ.npy test/hg38/ENCFF074VQD.npy \
		--chrom chr20 --out-file test/hg38/ENCFF622DXZ.ENCFF074VQD.score.txt
	```

4) Output looks like: (header: bootstrap_index, mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1, mseprom, msegene, mseenh).
	```bash
	bootstrap_-1	20.45688606636623	1730.3503548526915	195.52252657980728	0.01705378703206674	848	3462	2976	0.5852748736100822	0.590682173511888	376.1018309950674	31.24613030186926	94.01719916101615
	```


## Validation for submissions

In order to validate your BIGWIG. Use `validate.py`.

```bash
$ python validate.py [YOUR_SUBMISSION_BIGWIG]
```

## Ranking for submissions

1) [Generate bootstrap label](#how-to-generate-bootstrap-labels)

2) In order to speed up scoring, convert `TRUTH_BIGWIG` into numpy array/object (binned at `25`). Repeat this for each pair of cell type and assay.
	```bash
	$ python bw_to_npy.py [TRUTH_BIGWIG] --out-npy-prefix [TRUTH_NPY_PREFIX]
	```

3) Create a score database.
	```bash
	$ python db.py [SCORE_DB_FILE]
	```

4) For each assay type, build a variance `.npy` file, which calculates a variance for each bin for each chromosome across all cell types. Without this variance file, `msevar` will be `0.0`.
	$ python build_var_npy.py [TRUTH_NPY_CELL1] [TRUTH_NPY_CELL2] ... \
		--out-npy-prefix var_[ASSAY_OR_MARK_ID]

5) Score each submission with bootstrap labels. `--validated` is only for validated submissions binned at `25`. With this flag turned on, `score.py` will skip interpolation of intervals in a bigwig. For ranking, you need to define all metadata for a submission like `--cell [CELL_ID] --assay [ASSAY_OR_MARK_ID] -t [TEAM_ID_INT] -s [SUBMISSION_ID_INT]`. These values will be written to a database file together with bootstrap scores. Repeat this for each submission (one submission per team for each pair of cell type and assay).
	```bash
	$ python score.py [YOUR_VALIDATED_SUBMISSION_BIGWIG_OR_NPY] [TRUTH_NPY] \
	    --var-npy var_[ASSAY_OR_MARK_ID].npy \
		--db-file [SCORE_DB_FILE] \
		--cell [CELL_ID] --assay [ASSAY_OR_MARK_ID] \
		-t [TEAM_ID_INT] -s [SUBMISSION_ID_INT] \
		--validated
	```

6) Calculate ranks based on DB file
	```bash
	$ python rank.py [SCORE_DB_FILE]
	```


## For challenge admins

### How to generate bootstrap labels?

Download `submission_template.bigwig` from Synapse imputation challenge site. The following command will make 10-fold (default) bootstrap index for each chromosome. Output is a single `.npy` file which have all bootstrap labels for corresponding bootstrap index and chromosomes.

```bash
$ python build_bootstrapped_label.py submission_template.bigwig
```

### How to use bootstrapped label?

Simply run `score.py` with `--bootstrapped-label-npy bootstrapped_label.npy`.

