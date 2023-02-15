<img src="https://github.com/3DGenomes/metawaffle/blob/master/Logo_capture.png" height= "225" width="500">

This add-on of TADbit allows you to deconvolve the structural signal from a set of regions of interest and obtain the different structural clusters.  

<!-- TOC depthFrom:1 depthTo:8 withLinks:1 updateOnSave:1 orderedList:0 -->

  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Quick start](#quick-start)
      - [Preparing the input files](#preparing-the-input-files)
        - [BAM file](#bam-file)
        - [ChIP-peak file](#chip-peak-file)
      - [Extracting submatrices](#extracting-submatrices)
      - [Self-organizing feature map](#self-organizing-feature-map)
  - [Usage](#usage)
  - [Contributors](#contributors)
  - [Citation](#citation)

<!-- /TOC -->

# Requirements
Metawaffle was written and tested using Python 2.7.
It requires the packages `TADbit` (0.4.97),`Neupy` (0.6.4),`matplotlib` (3.2.1), `scipy` (1.0.0), `numpy` (1.14.0), `pandas` (0.22.0),`pysam` (0.15.2), `pickle` (2.0.0), `functools` (3.2.3), `multiprocess` (0.70.9) and `scikit-learn` (0.22.2).

# Installation

In order to install the package you can download the code from our [GitHub repo](https://github.com/3DGenomes/metawaffle) and install it manually. If needed, the dependencies will be downloaded automatically; if you encounter problems, try to install the other [requirements](#requirements) manually.

```bash
git clone https://github.com/3DGenomes/metawaffle # or download manually
cd metawaffle
python setup.py install
```

# Quick start
After the installation, you can run the provided example to familiarize with the functions of metawaffle.

Use `metawaffle -h` for quick help and orientation.
## Preparing the input files

### BAM file
For the tutorial will be used the Hi-C data from the chromosome 22 of GM12878 cell line from the merged replicates listed in (META-WAFFLE ARTICLE).

The biases have to be computed before running `bam2count`. We use `TADbit` as it gives the possibility to obtain the biases using various normalization methods (check https://github.com/3DGenomes/TADbit). For this example we have extracted the biases using the OneD normalization (E Vidal et al. 2018) at 5 kb resolution. The biases will be stored in a pickle file. We will use a 5 kb resolution.

```bash
metawaffle bam2count \
example/chr22_GM12878.bam \
5000 \
example/ \
example/biases_5kb_GM12878.pickle \
```
You can notice that if you do not have a BAM format file, you can use a tabseparated format with the same information as the output file. The output file contains 4 columns: 

- bin i: genomic bin i.
- bin j: genomic bin j.
- raw contact value: The raw interaction value.
- normalized contact value: The normalized interaction value.

### ChIP-peak file
For the tutorial we will use the ChIP-seq peaks from CTCF in the chromosome 2 (listed in PAPER).First of all we can run a quality check to know the distribution and the width of the peaks.

```bash
metawaffle check_peaks \
example/ctcf_chr2.bed \
example/ \
```
The quality check will provide worthy information to decide the genomic interval of distance and the windows span around the center of the peak. Then, the pair list will be generated running the `pairlist` script.

<img src="https://github.com/3DGenomes/metawaffle/blob/master/images/test_evaluation.png" height= "350" width="470">

```bash
metawaffle pairlist \
example/ctcf_chr2.bed \
example/ \
example/human_grch38_size \
ctcf_chr2 \
27500 \
2500000 \
255000-2500000 \
```
This command will generate a file with the coordinate pairs with the CTCF peak centered in 55 kb region. Only those coordinates within the interval distance of 255 kb - 2.5 Mb will be stored.

The output will be a two column file:

*chr22:40502977-40547977	chr22:41644314-41689314*

*chr22:39308838-39353838	chr22:40502977-40547977*

*chr22:40502977-40547977	chr22:41650826-41695826*

*chr22:39291974-39336974	chr22:40502977-40547977*

*chr22:39492032-39537032	chr22:40502977-40547977*


## Extracting submatrices
The next step will be to extract the contact matrices from the coordinate pair list. For the example we will extract the all the matrices using the tag -A True. If you only want to use a random set of coordinates from the pairlist, -A True and specify the sample size.

```bash
metawaffle peak2matrix \
example/ctcf_chr2_255000_2500000.tsv \
example/human_grch38_size \
5000 \
example/tmpdir/ \
example/ \
8 \
ctcf_chr2 \
example/biases_5kb_GM12878.pickle \
example/matrix/ \
27500 \
-A True \
```
The output file will contain row-wise the coordinates and their contact matrices from the pairlist file.

## Self-organizing feature map
To classify the extracted contact matrices according to their structural pattern, a competitive neural network provided in `Neupy` package, the Self-Organizing Feature Map (SOFM or SOM).

```bash
metawaffle sofm \
example/matrices_ctcf_chr2.tsv \
example/sofm/ \
ctcf_chr2 \
20 \
12 \
5 \
1 \
0.01 \
0.01 \
-plot True \
```
If you want to know more about SOFM: http://neupy.com/apidocs/neupy.algorithms.competitive.sofm.html .
After running this command, you will obtain various files:

- Cluster_n.bed: You will obtain bed files for each of the neurons according to the grid size, with the classified coordinates.
- Heatmap_info : image with the number of coordinates classified per neuron, in order to check how the coordinates have been distributed in the SOFM map.
<img src="https://github.com/3DGenomes/metawaffle/blob/master/images/heatmap_info_ctcf_chr2.png" height= "350" width="350">
- weight_matrix.npy: In here you will have the weight of each of the input values provided to SOFM.
- sofm.pickle: The model trained. In order to repeat or re-run the analysis.
- If -plot True: image of the SOFM map. 
<img src="https://github.com/3DGenomes/metawaffle/blob/master/images/ctcf_chr2_SOFM.png" height= "350" width="350">


# Usage
metawaffle has these basic commands: 
* `bam2count` to convert the bamfile into tabseparated counts.
* `check_peaks` to evaluate peaks of interest, and get the widht and distance between contiguous peaks.
* `pairlist` used to generate the coordinates list to extract the matrices.
* `peak2matrix` contains the tools to extract the contact matrices from the pairs list.
* `sofm` used to classify the contact matrices according to their structural pattern.


Frequently asked questions
--------------------------

Check the label `FAQ <https://github.com/3DGenomes/TADbit/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3AFAQ+>`_ in TADbit issues.

If your question is still unanswered feel free to open a new issue.

# Contributors
This add-on of TADbit is currently developed at the  `MarciusLab <http://www.marciuslab.org>`_ with the contributions of Silvia Galan and François Serra.

# Citation
Please, cite this article if you use TADbit.

Serra, F., Baù, D., Goodstadt, M., Castillo, D. Filion, G., & Marti-Renom, M.A. (2017).
**Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors.**
*PLOS Comp Bio* 13(7) e1005665. `doi:10.1371/journal.pcbi.1005665 <https://doi.org/10.1371/journal.pcbi.1005665>`_

My successfull run
--------------------------
Some notes before run: 
- temp directory must be empty before run! (otherwise it returns *ValueError: too many values to unpack*)
- pickle file for biases was generated by TADbit's oneD function

pairlist function: 
metawaffle pairlist ../data/CutNTag/bed/promoters_CTCF_no_G4_raw.bed output/ ../data/mm10.chrom.sizes.txt prom_ctcf_chr1 25000 2500000 255000-2500000

Notes:
- It creates CTCF-CTCF cis pairs
- Extraction of Hi-C submatrices of all possible cis pairs of CTCF peaks linearly separated by between 255 kb and 2.5 Mb
- Submatrix dimension: one range of the pair (row) spans 25 kb (10 kb upstream and 10 kb downstream with the central 5 kb bin containing a CTCF peak)

peak2matrix function:
metawaffle peak2matrix output/ctcf_chr1_255000_2500000.tsv ../data/mm10.chrom.sizes.txt 10000 output/tmpdir/ output/ 8 ctcf_chr1 TADbit_normalization/biases_10kb.pickle output/matrix/ 25000 -A True

Notes: 
- output/matrix/ folder must contain the HiC interaction matrix with 10 kb resolution 
	Here I used a **HiC matrix divided by 10 kb** (that was the resolution I used for the peak2matrix function!). 
	In addition, I made two helper functions in the peak2matrix function for the run (*within_region, around_region* functions).
	*around_region* function enlarges the HiC regions with 20 for both sides.  
- TADbit_normalization/biases_10kb.pickle - this pickle object comes from the oneD normalization made by TADbit
- I made two helper functions in the peak2matrix function (peak2matrix .py file).

SOFM function: 
metawaffle sofm output/matrices_ctcf_chr1.tsv output/sofm/ ctcf_chr1 30 6 30 5 0.5 0.01 -plot True
