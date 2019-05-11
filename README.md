# ENCODE Imputation Challenge Scoring and Validation Scripts

## Installation for dependencies (pip)

```bash
$ pip install numpy scikit-learn pyBigWig
```

## Installation for dependencies (Conda)

1) [Install Conda](https://conda.io/docs/user-guide/install/linux.html) first.

2) Install `numpy`, `scikit-learn` and `pyBigWig`.
	```bash
	$ conda install -c bioconda numpy scikit-learn pyBigWig 
	```

## Example (hg38)

1) Download ENCFF622DXZ and ENCFF074VQD from ENCODE portal.
	```bash
	$ mkdir -p test/hg38 && cd test/hg38
	$ wget https://www.encodeproject.org/files/ENCFF622DXZ/@@download/ENCFF622DXZ.bigWig
	$ wget https://www.encodeproject.org/files/ENCFF074VQD/@@download/ENCFF074VQD.bigWig
	```

2) Run it.
	```bash
	$ python score.py test/hg38/ENCFF622DXZ.bigWig test/hg38/ENCFF074VQD.bigWig \
		--gene-annotations annot/hg38/gencode.v29.genes.gtf.bed.gz \
		--enh-annotations annot/hg38/F5.hg38.enhancers.bed.gz \
		--chrom chr20 \
		--window-size 25 --prom-loc 80 --nth 1 \
		--out test/hg38/ENCFF622DXZ.ENCFF074VQD.score.txt
	```

3) Output looks like: (header: chrom, mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1, mseprom, msegene, mseenh).
	```bash
	chr20	20.45688606636623	1730.3503548526915	195.52252657980728	0.01705378703206674	848	3462	2976	0.5852748736100822	0.590682173511888	376.1018309950674	31.24613030186926	94.01719916101615
	```


## Validation for submissions

In order to validate your BIGWIG. Use `validate.py`.

```bash
$ python validate.py [YOUR_SUBMISSION_BIGWIG]
```
