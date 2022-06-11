# Pynoma

![](https://img.shields.io/badge/python-v3.x-blue)

## Summary

- [Introduction](#introduction)
- [Installation](#installation)
- [Search Types](#search-types)
    - [Search by gene](#search-by-gene)
    - [Search by region](#search-by-region)
    - [Search by transcript](#search-by-transcript)
    - [Search by variant](#search-by-variant)
- [Batch search](#batch-search)
- [BibTeX entry](#bibtex-entry) 
- [Acknowledgement](#acknowledgement)

## Introduction

Pynoma is an API developed to facilitate the access to human variant data from gnomAD database, working both with gnomAD version 2 and 3.
The package retrieves both regular as well as clinical data, and offers support to four kinds of different searches, as well as the possibility to search in batches.
For plotting the data, please take a look at the [BIOVARS package](https://github.com/bioinfo-hcpa/biovars). If you have scientific interests or want to use our package in formal reports, we kindly ask you to cite us in your publication: [Carneiro, P., Colombelli, F., Recamonde-Mendoza, M., and Matte, U. (2022). Pynoma, PyABraOM and BIOVARS: Towards genetic variant data acquisition and integration. bioRxiv.](#bibtex-entry)

## Installation


Currently there is not a PyPI version for the Pynoma API, so the installation needs that you clone this repository and install it as local package.

    $ git clone https://github.com/bioinfo-hcpa/pynoma.git
    $ pip install -e pynoma

## Search Types

There are 4 kinds of different searches supported by pynoma. Gene, transcript and region searches return the same kind of output, while variant searches build a different dataframe format.

Since gnomAD have two versions, the user must specify which one is to be used by choosing an integer value of either 2 or 3. Additionally, when searching for genes, transcripts or a region in chromosomes X or Y, a new column "Number of Hemizygotes" will be added to the outputted dataframe, so the user should have caution when performing pandas concatenation operations or batch searchings that could potentially mix both kinds of dataframe, resulting in table cells with NaN values.

### Search by gene

GeneSearch(gnomad_version: int, gene: str)<br />
.get_data(standard=True, additional_population_info=False)

```python
from pynoma import GeneSearch
gs = GeneSearch(3, "IDUA")
df, clinical_df = gs.get_data()
```

### Search by transcript


TranscriptSearch(gnomad_version: int, transcript_id: str)<br />
.get_data(standard=True, additional_population_info=False)

```python
from pynoma import TranscriptSearch
ts = TranscriptSearch(3, "ENST00000247933")
df, clinical_df = ts.get_data()
```


### Search by region

RegionSearch(gnomad_version: int, chromosome, region_start: int, region_end: int)<br />
.get_data(standard=True, additional_population_info=False)

```python
from pynoma import RegionSearch
rs = RegionSearch(3, 4, 1002741, 1002771)
df, clinical_df = rs.get_data()
```

### Search by variant

VariantSearch(gnomad_version: int, variant_id: str)<br />
.get_data(raw=False)

```python
from pynoma import VariantSearch
vs = VariantSearch(3, '4-1002747-G-A')
df, meta = vs.get_data()
```

## Batch search

If the user wants to configure multiple searches, including different ones (gene, transcript, region) with the exception of variant searches (that have different dataframe formats), they can use the batch search function.
For example, to search for variants in a list of 5 genes, let's say _ACE2_, _BRCA_, _ID4_, _MTOR_ and _EMP1_, the batch search can be used as follows:

```python
from pynoma import helper, GeneSearch
genes = [GeneSearch(3, "ACE2"), GeneSearch(3, "BRCA"), GeneSearch(3, "ID4"), GeneSearch(3, "MTOR"), GeneSearch(3, "EMP1")]
df = helper.batch_search(genes, standard=True, additional_population_info=False, verbose=True)
```

Besides the list of Search objects, the other parameters (standard, additional_population_info and verbose) follow the same logic of the individual searches.


## BibTeX entry

```
@article {Carneiro2022.06.07.495190,
	author = {Carneiro, Paola and Colombelli, Felipe and Recamonde-Mendoza, Mariana and Matte, Ursula},
	title = {Pynoma, PyABraOM and BIOVARS: Towards genetic variant data acquisition and integration},
	elocation-id = {2022.06.07.495190},
	year = {2022},
	doi = {10.1101/2022.06.07.495190},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.07.495190},
	eprint = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.07.495190.full.pdf},
	journal = {bioRxiv}
}
```

## Acknowledgement

This research was supported by the National Council for Scientific and Technological Development (CNPq) and the Research Incentive Fund (FIPE) from Hospital de Clínicas de Porto Alegre.
