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

## Example (hg19)

1) Download ENCFF141LFV and ENCFF507RBF from ENCODE portal.
	```bash
	$ mkdir -p test/hg19 && cd test/hg19
	$ wget https://www.encodeproject.org/files/ENCFF141LFV/@@download/ENCFF141LFV.bigWig
	$ wget https://www.encodeproject.org/files/ENCFF507RBF/@@download/ENCFF507RBF.bigWig
	```

2) Run it.
	```bash
	$ python score.py test/hg19/ENCFF141LFV.bigWig test/hg19/ENCFF507RBF.bigWig \
		--gene-annotations annot/hg19/gencode.v19.annotation.protein_coding.full.sorted.genes.bed.gz \
		--enh-annotations annot/hg19/human_permissive_enhancers_phase_1_and_2.bed.gz \
		--chrom chr21 chr22 \
		--window-size 25 --prom-loc 80 --nth 1 \
		--out test/hg19/ENCFF141LFV.ENCFF507RBF.score.txt
	```

3) Output looks like: (header: chrom, mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1, mseprom, msegene, mseenh).
	```bash
	$ cat ENCFF141LFV.ENCFF507RBF.score.txt
	chr21	0.7319117871666568	3.2720607489506124	37.598248674801255	0.2399089043479622	1188	3894	3266	0.7485369355214815	0.6948649403632234	3.236683554034965	1.421404590155849	6.032906435152852
	chr22	1.366929166612142	5.13244903384441	61.82154443671531	0.29302792900228364	1436	5645	4589	0.8086374687050369	0.7275331625264598	4.609551163010887	2.3647588930301713	7.901332816276359
	```
	**chrom**|**mse**|**mse1obs**|**mse1imp**|**gwcorr**|**match1**|**catch1obs**|**catch1imp**|**aucobs1**|**aucimp1**|**mseprom**|**msegene**|**mseenh**
	:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
	chr21|0.7319117871666568|3.2720607489506124|37.598248674801255|0.2399089043479622|1188|3894|3266|0.7485369355214815|0.6948649403632234|3.236683554034965|1.421404590155849|6.032906435152852
	 |chr22|1.366929166612142|5.13244903384441|61.82154443671531|0.29302792900228364|1436|5645|4589|0.8086374687050369|0.7275331625264598|4.609551163010887|2.3647588930301713|7.901332816276359
