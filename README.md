# imputation_challenge
ENCODE Imputation Challenge scoring &amp; validation scripts

## Installation (pip)

```bash
$ pip install numpy scikit-learn pyBigWig
```

## Installation (Conda)

1) [Install Conda](https://conda.io/docs/user-guide/install/linux.html) first.

2) Install `numpy`, `scikit-learn` and `pyBigWig`.
	```bash
	$ conda install -c bioconda numpy scikit-learn pyBigWig 
	```

## Example

```bash
$ python score.py ENCFF867EID.bigWig ENCFF053MKO.bigWig \
	--gene-annotations gencode.v19.annotation.protein_coding.full.sorted.genes.bed \
	--enh-annotations human_permissive_enhancers_phase_1_and_2.bed \
	--chrom chr3 chr21 \
	--window-size 25 --prom-loc 80 --nth 1 \
	--out ENCFF053MKO.ENCFF867EID.score.txt
```