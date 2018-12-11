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

## Example (hg19)

1) Download ENCFF141LFV and ENCFF507RBF from ENCODE portal.
	```bash
	$ wget https://www.encodeproject.org/files/ENCFF141LFV/@@download/ENCFF141LFV.bigWig
	$ wget https://www.encodeproject.org/files/ENCFF507RBF/@@download/ENCFF507RBF.bigWig
	```

2) Download `gencode.v19.annotation.protein_coding.full.sorted.genes.bed` and `human_permissive_enhancers_phase_1_and_2.bed`.

3) Run it.
	```bash
	$ python score.py ENCFF141LFV.bigWig ENCFF507RBF.bigWig \
		--gene-annotations gencode.v19.annotation.protein_coding.full.sorted.genes.bed \
		--enh-annotations human_permissive_enhancers_phase_1_and_2.bed \
		--chrom chr21 chr22 \
		--window-size 25 --prom-loc 80 --nth 1 \
		--out ENCFF141LFV.ENCFF507RBF.score.txt
	```

4) Output looks like: (header: chrom, mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1, mseprom, msegene, mseenh).
	```bash
	$ cat 
	chr21	0.7319117871666568	3.2720607489506124	37.598248674801255	0.2399089043479622	1188	3894	3266	0.7485369355214815	0.6948649403632234	3.236683554034965	1.421404590155849	6.032906435152852
	chr22	1.366929166612142	5.13244903384441	61.82154443671531	0.29302792900228364	1436	5645	4589	0.8086374687050369	0.7275331625264598	4.609551163010887	2.3647588930301713	7.901332816276359
	```
	**chrom**|**mse**|**mse1obs**|**mse1imp**|**gwcorr**|**match1**|**catch1obs**|**catch1imp**|**aucobs1**|**aucimp1**|**mseprom**|**msegene**|**mseenh**
	:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
	chr21|0.7319117871666568|3.2720607489506124|37.598248674801255|0.2399089043479622|1188|3894|3266|0.7485369355214815|0.6948649403632234|3.236683554034965|1.421404590155849|6.032906435152852
	 |chr22|1.366929166612142|5.13244903384441|61.82154443671531|0.29302792900228364|1436|5645|4589|0.8086374687050369|0.7275331625264598|4.609551163010887|2.3647588930301713