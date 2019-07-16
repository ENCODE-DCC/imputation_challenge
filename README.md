# ENCODE Imputation Challenge Scoring/Ranking and Validation Scripts

## Installation

1) [Install Conda](https://docs.conda.io/en/latest/miniconda.html) first.

2) Install `numpy`, `scikit-learn` and `pyBigWig`.
	```bash
	$ conda install -y -c bioconda numpy scikit-learn pyBigWig sqlite scipy
	```

## Validating a submission

```bash
$ python validate.py [YOUR_SUBMISSION_BIGWIG]
```

## Scoring a submission

1) Download ENCFF622DXZ and ENCFF074VQD from ENCODE portal.
	```bash
	$ mkdir -p test/hg38 && cd test/hg38
	$ wget https://www.encodeproject.org/files/ENCFF622DXZ/@@download/ENCFF622DXZ.bigWig
	$ wget https://www.encodeproject.org/files/ENCFF074VQD/@@download/ENCFF074VQD.bigWig
	```

2) Convert it to numpy array. This is to speed up scoring multiple submissions. `score.py` can also take bigwigs so you can skip this step.
	```bash
	$ python bw_to_npy.py test/hg38/ENCFF622DXZ.bigWig
	$ python bw_to_npy.py test/hg38/ENCFF074VQD.bigWig
	```

3) Run it. If you score without a variance `.npy` file specified as `--var-npy`, then `msevar` metric will be `0.0`.
	```bash
	$ python score.py test/hg38/ENCFF622DXZ.npy test/hg38/ENCFF074VQD.npy --chrom chr20
	```

## Ranking for submissions

1) Create a score database.
	```bash
	$ python db.py [NEW_SCORE_DB_FILE]
	```

2) In order to speed up scoring, convert `TRUTH_BIGWIG` into numpy array/object (binned at `25`). Repeat this for each pair of cell type and assay. `--out-npy-prefix [TRUTH_NPY_PREFIX]` is optional. Repeat this for all truth bigwigs.
	```bash
	$ python bw_to_npy.py [TRUTH_BIGWIG] --out-npy-prefix [TRUTH_NPY_PREFIX]
	```

3) For each assay type, build a variance `.npy` file, which calculates a variance for each bin for each chromosome across all cell types. Without this variance file, `msevar` will be `0.0`.
	```bash
	$ python build_var_npy.py [TRUTH_NPY_CELL1] [TRUTH_NPY_CELL2] ... --out-npy-prefix var_[ASSAY_OR_MARK_ID]
	```

4) Score each submission. `--validated` is only for a validated bigwig submission binned at `25`. With this flag turned on, `score.py` will skip interpolation of intervals in a bigwig. For ranking, you need to define metadata for a submission like -t [TEAM_ID_INT] -s [SUBMISSION_ID_INT]`. These values will be written to a database file together with bootstrap scores. Repeat this for each submission (one submission per team for each pair of cell type and assay).
	```bash
	$ python score.py [YOUR_VALIDATED_SUBMISSION_BIGWIG_OR_NPY] [TRUTH_NPY] \
	    --var-npy var_[ASSAY_OR_MARK_ID].npy \
		--db-file [SCORE_DB_FILE] \
		--validated \
		-t [TEAM_ID_INT] -s [SUBMISSION_ID_INT]
	```

5) Calculate ranks based on DB file
	```bash
	$ python rank.py [SCORE_DB_FILE]
	```

## Setting up a leaderboard server (admins only)

1) Create a server instance on AWS.

2) Install Synapse client.
	```bash
	$ pip install synapseclient
	```

3) Authenticate yourself on the server
	```bash
	$ synapse login --remember-me -u [USERNAME] -p [PASSWORD]
	```

4) Create a score database.
	```bash
	$ python db.py [NEW_SCORE_DB_FILE]
	```

5) Run `score_leaderboard.py` for each submission. Files on `TRUTH_NPY_DIR` should looks like `CXXMYY.npy`. Files on `VAR_NPY_DIR` should look like `var_MXX.npy`. Submissions will be downloaded on `SUBMISSION_DOWNLOAD_DIR`.
	```bash
	$ python score_leaderboard.py [EVALUATION_QUEUE_ID] [TRUTH_NPY_DIR] \
	    --var-npy-dir [VAR_NPY_DIR] \
	    --submission-dir [SUBMISSION_DOWNLOAD_DIR] \
	    --send-to-admin \
		--db-file [SCORE_DB_FILE]
	```

	Example:
	```bash
	$ python score_leaderboard.py 9614278 --send-to-admin --submission-dir ./submissions --db-file test.db
	```
