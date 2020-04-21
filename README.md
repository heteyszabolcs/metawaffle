<img src="https://github.com/3DGenomes/meta-waffle/blob/master/Logo_capture.png" height= "225" width="500">

This add-on of TADbit allows you to deconvolve the structural signal from a set of regions of interest and obtain the different structural clusters.  

meta-waffle is a tool for the extraction and classification of 

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
Meta-waffle was written and tested using Python 3.6, 3.7.
It requires the packages `TADbit` (0.4.97),`Neupy` (0.6.4),`matplotlib` (3.2.1), `scipy` (1.0.0), `numpy` (1.14.0), `pandas` (0.22.0),`pysam` (0.15.2), `pickle` (2.0.0), `functools` (3.2.3), `multiprocess` (0.70.9) and `scikit-learn` (0.22.2).

# Installation

In order to install the package you can download the code from our [GitHub repo](https://github.com/3DGenomes/meta-waffle) and install it manually. If needed, the dependencies will be downloaded automatically; if you encounter problems, try to install the other [requirements](#requirements) manually.

```bash
git clone https://github.com/3DGenomes/meta-waffle # or download manually
cd meta-waffle
python setup.py install
```

# Quick start
After the installation, you can run the provided example to familiarize with the functions of meta-waffle.

Use `metawaffle -h` for quick help and orientation.
## Preparing the input files

### BAM file
For the tutorial will be used the Hi-C data from the chromosome 22 of GM12878 cell line from the merged replicates listed in (META-WAFFLE ARTICLE).

The biases have to be computed before running `bam2count`. We use `TADbit` as it gives the possibility to obtain the biases using various normalization methods (check https://github.com/3DGenomes/TADbit). For this example we have extracted the biases using the OneD normalization (E Vidal et al. 2018) at 5 kb resolution. The biases will be stored in a pickle file. We will use a 5 kb resolution.

```bash
metawaffle bam2count \
-bam examples/chr22_GM12878.bam \
-r 5000 \
-o examples/ \
-b examples/biases_chr22_GM12878.pickle \
```
You can notice that if you do not have a BAM format file, you can use a tabseparated format with the same information as the output file. The output file contains 4 columns: 

- bin i: genomic bin i.
- bin j: genomic bin j.
- raw contact value: The raw interaction value.
- normalized contact value: The normalized interaction value.

### ChIP-peak file
For the tutorial we will use the ChIP-seq peaks from CTCF in the chromosome 22 (listed in PAPER).First of all we can run a quality check to know the distribution and the width of the peaks.

```bash
metawaffle check_peaks \
-i CTCF_chr22.bed \
-o examples/ \
```
The quality check will provide worthy information to decide the genomic interval of distance and the windows span around the center of the peak. Then, the pair list will be generated running the `pairlist` script.

![Image](https://github.com/3DGenomes/meta-waffle/tree/master/images/test_evaluation.png)


```bash
metawaffle pairlist \
-i examples/CTCF_chr22.bed \
-b examples/chr22_GM12878.bam \
-o examples/ \
-n ctcf_example \
-s  9 \
-m 1000000 \
-w 255000-1000000 \
```
[OUTPUT!!]

## Extracting submatrices
The next step will be to extract the contact matrices from the coordinate pair list.
```bash
metawaffle peak2matrix \
-i examples/CTCF_chr22_pairs.bed \
-bam examples/chr22_GM12878.bam \
-r 5000 \
-t examples/tmpdir/ \
-o examples/ \
-C 8 \
-n ctcf_example \
-b examples/biases_chr22_GM12878.pickle \
-mats 
-m 1000000 \
-w 255000-1000000 \
-specific True \
```

[OUTPUT!!]

## Self-organizing feature map
To classify the extracted contact matrices according to their structural pattern, a competitive neural network provided in `Neupy` package, the Self-Organizing Feature Map (SOFM or SOM).

```bash
metawaffle sofm \
-i examples/CTCF_chr22_pairs.bed \
-o examples/ \
-n ctcf_example \
-sigma 0.3 \
-e 10 \
-s 9 \
-g 4 \
-l 0.2 \
-std \
-step \
```
[OUTPUT!!]

# Usage
meta-waffle has these basic commands: 
* `bam2count` to convert the bamfile into tabseparated counts.
* `pairlist` used to generate the coordinates list to extract the matrices.
* `p2m` contains the tools to extract the contact matrices from the pairs list.
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
