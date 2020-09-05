in_region = """query VariantInRegion($chrom: String!, $start: Int!, $stop: Int!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
  region(start: $start, stop: $stop, chrom: $chrom, reference_genome: $referenceGenome) {
    clinvar_variants {
      clinical_significance
      clinvar_variation_id
      gold_stars
      major_consequence
      pos
      variant_id
    }
    variants(dataset: $datasetId) {
      consequence
      flags
      gene_id
      gene_symbol
      hgvs
      hgvsc
      hgvsp
      lof
      lof_filter
      lof_flags
      pos
      rsid
      variant_id: variantId
      exome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
      genome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
    }
  }
}"""

in_region_variables = """
            {
                "chrom":"%s",
                "datasetId":"%s",
                "referenceGenome":"GRCh38",
                "start":%s,
                "stop":%s
            }"""
    
region_coverage = """query RegionCoverage($chrom: String!, $start: Int!, $stop: Int!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!, $includeExomeCoverage: Boolean!, $includeGenomeCoverage: Boolean!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referenceGenome) {
    exome_coverage(dataset: $datasetId) @include(if: $includeExomeCoverage) {
      pos
      mean
      median
      over_1
      over_5
      over_10
      over_15
      over_20
      over_25
      over_30
      over_50
      over_100
    }
    genome_coverage(dataset: $datasetId) @include(if: $includeGenomeCoverage) {
      pos
      mean
      median
      over_1
      over_5
      over_10
      over_15
      over_20
      over_25
      over_30
      over_50
      over_100
    }
  }
}"""
    
fetch_region = """query FetchRegion($chrom: String!, $start: Int!, $stop: Int!, $referenceGenome: ReferenceGenomeId!) {
    region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referenceGenome) {
      reference_genome
      chrom
      start
      stop
      genes {
        gene_id
        symbol
        start
        stop
        exons {
          feature_type
          start
          stop
        }
      }
    }
  }"""



gene_id = """query Search($dataset: DatasetId!, $query: String!) {
    searchResults(dataset: $dataset, query: $query) {
      label
      value: url
    }
  }"""


gene_id_variables = """{
  "dataset": "%s",
  "query": "%s"
}"""