# pynoma

![](https://img.shields.io/badge/python-v3.x-blue)

## Introduction


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

## Features

Not available yet.
A class Helper will be provided to attend:
* Batch searching 
* Data filtering methods
* Saving dataframe methods

## Additional Information

## References
