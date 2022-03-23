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

## Introduction

Pynoma is an API developed to facilitate the access to human variant data from gnomAD database, working both with gnomAD version 2 and 3.
The package retrieves both regular as well as clinical data, and offers support to four kinds of different searches, as well as the possibility to search in batches.
For plotting the data, please take a look at the [BIOVARS package](https://github.com/bioinfo-hcpa/biovars).

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
gs = GeneSearch(3, "idua")
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

