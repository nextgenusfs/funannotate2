[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate2.svg)](https://github.com/nextgenusfs/funannotate2/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate2)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# funannotate2: eukaryotic genome annotation pipeline


#### This is work in progress. Do not expect it to work until a release is tagged.

Setting up a conda environment with these packages is necessary, then it can be installed with pip.

```shell
mamba create -n funannotate2 "python<=3.10" biopython "evidencemodeler>=2" minimap2 miniprot snap "augustus==3.5.0" glimmerhmm diamond trnascan-se table2asn
conda activate funannotate2
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git
```

Additional tools like genemarkHMM must be installed manually due to licensing.

To install release versions use the pip package manager, like so:
```
python -m pip install funannotate2
```

To install the most up to date code from this repo, you can run:
```
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git --upgrade --force --no-deps
```
