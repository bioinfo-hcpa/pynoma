"""This module contains the Search class, which is used to search the gnomAD database."""
from typing import Union
from requests import post
import pandas as pd
from pynoma.DataManager import DataManager
from pynoma.Logger import Logger

class Search:


    dataset_id_map = {
        2: "gnomad_r2_1",
        "2": "gnomad_r2_1",
        "hg19": "gnomad_r2_1",
        "hg37": "gnomad_r2_1",
        3: "gnomad_r3",
        "3": "gnomad_r3",
        "hg38": "gnomad_r3"
    }

    def __init__(self, dataset_version: Union[int, str], query: str, query_variables: str):
        """Constructor for the Search class.

        Args:
            dataset_version: The version of the gnomAD dataset to be used. It can be either 2, or hg19/h38
            query: The query to be used to search the gnomAD database.
            query_variables: The variables to be used in the query.
        """
        self.end_point = "https://gnomad.broadinstitute.org/api/"
        self.query = query
        self.query_vars = query_variables
        self.dataset_id, self.reference_genome = self.get_dataset_id(dataset_version)

        self.dm = None   # attribute holding DataManager object

    
    def request_gnomad(self, variables: Union[str, tuple]) -> dict:
        """Send a POST request to the gnomAD API.

        Args:
            variables: The variables to be used in the query. See examples in the Queries file.

        Returns:
            The response JSON from the gnomAD API request.
        """
        variables = self.query_vars % variables
        response = post(self.end_point, data={'query': self.query, 'variables': variables}, timeout=None)
        return response.json()

    
    @classmethod
    def get_dataset_id(cls, version: Union[int, str]):
        """Get the dataset ID based on the user version provided.

        The user can provide the gnomAD version as an integer (2 or 3) or the reference genome as a string ("hg19" or
        "hg38"). For now, the only available versions are gnomAD's version 2 and 3.

        Args:
            version: The version of the gnomAD dataset to be used. It can be either 2, 3 or hg19/h38

        Raises:
            Exception: If the version provided is not valid.
        
        Returns:
            The dataset ID and the reference genome.
        """
        if version in cls.dataset_id_map:
            dataset_id = cls.dataset_id_map[version]
            if dataset_id == "gnomad_r2_1":
                reference_genome = "GRCh37"
            else:
                reference_genome = "GRCh38"
        else:
            raise Exception("Invalid dataset version. Choose gnomAD's version 2, 3 or hg19/h38.")
        return dataset_id, reference_genome


class RegionSearch(Search):

    def __init__(self, 
                 dataset_version: Union[int, str],
                 chromosome: Union[int, str], 
                 start_position: Union[int, str], 
                 end_position: Union[int, str]):
        """Constructor for the RegionSearch class.
        
        Args:
            dataset_version: The version of the gnomAD dataset to be used. It can be either 2, 3 or hg19/h38
            chromosome: The chromosome number to search for.
            start_position: The start position of the region to search for.
            end_position: The end position of the region to search for.
        """
        from pynoma.Queries import in_region_v3, in_region_v2, in_region_variables
        super().__init__(dataset_version, "", in_region_variables)
        dataset_id, _ = self.get_dataset_id(dataset_version)
        if dataset_id == "gnomad_r2_1":
            self.query = in_region_v2
        else:
            self.query = in_region_v3

        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)


    def get_json(self) -> dict:
        """Get the JSON data from the gnomAD API."""
        variables = (self.chromosome, self.dataset_id, self.reference_genome, self.start, self.end)
        return self.request_gnomad(variables)

    
    def get_data(self, 
                 standard=True, 
                 additional_population_info=False
                 ) -> tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
        """Get the data from the gnomAD API.

        Args:
            standard: If True, the data will be processed and returned in a standard format. Defaults to True.
            additional_population_info: If True, 9 additional columns with population information will be added to the
                dataframe. This information includes the allele frequency for each variant in 9 different populations.
                Defaults to False.

        Returns:
            A tuple containing the dataframe with the data and the clinical dataframe.
        """
        json_data = self.get_json()
        if not json_data['data']['region']['variants']:
            Logger.no_variants_found()
            return (None, None)

        self.dm = DataManager(json_data, self.dataset_id)

        if standard:
            self.dm.process_standard_dataframe()
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('standard'), self.dm.clinical_df
            return self.dm.standard_df, self.dm.clinical_df  # TODO: investigate linter type-checking complaint
        
        else:
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('raw'), self.dm.clinical_df
            return self.dm.raw_df, self.dm.clinical_df




class GeneSearch(Search):

    def __init__(self, dataset_version: Union[int, str], gene: str):

        from pynoma.Queries import gene_id, gene_id_variables
        super().__init__(dataset_version, gene_id, gene_id_variables)
        
        self.gene = gene
        self.gene_ens_id = None
        if not self.get_ensembl_id(gene_id_variables):
            return 

        from pynoma.Queries import variant_in_gene, variant_in_gene_variables
        self.query = variant_in_gene
        self.query_vars = variant_in_gene_variables


    def get_ensembl_id(self, query_variables):
        json_data = self.request_gnomad((self.gene, self.reference_genome))
        if not json_data['data']['gene_search']:
            Logger.no_gene_found_with_given_name(self.gene)
            return False
        
        self.gene_ens_id = json_data['data']['gene_search'][0]['ensembl_id']
        return True


    def get_gene_information(self):
        gene_info = self.request_gnomad(self.gene_ens_id)
        self.chromosome = gene_info['data']['gene']['chrom']
        self.start = gene_info['data']['gene']['start']
        self.end = gene_info['data']['gene']['stop']
        return

    def get_json(self):
        variables = (self.dataset_id, self.gene_ens_id)
        return self.request_gnomad(variables)


    def get_data(self, standard=True, additional_population_info=False):
        if not self.gene_ens_id:
            return (None, None)
        json_data = self.get_json()
        if not json_data['data']['gene']['variants']:
            Logger.no_variants_found()
            return (None, None)

        self.dm = DataManager(json_data, self.dataset_id, second_level_key='gene')

        if standard:
            self.dm.process_standard_dataframe()
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('standard'), self.dm.clinical_df
            return self.dm.standard_df, self.dm.clinical_df
        
        else:
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('raw'), self.dm.clinical_df
            return self.dm.raw_df, self.dm.clinical_df



class TranscriptSearch(Search):
    def __init__(self, dataset_version:int, transcript: str):
        from pynoma.Queries import variant_in_transcript, variant_in_transcript_variables
        
        self.transcript = transcript
        super().__init__(dataset_version, variant_in_transcript, variant_in_transcript_variables)
        
        
    def get_data(self, standard=True, additional_population_info=False):
        json_data = self.get_json()

        if not json_data['data']['transcript']['variants']:
            Logger.no_variants_found_for_given_transcript(self.transcript)
            return (None, None)

        json_data['data']['region'] = json_data['data'].pop('transcript')
        self.dm = DataManager(json_data, self.dataset_id)

        if standard:
            self.dm.process_standard_dataframe()
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('standard'), self.dm.clinical_df
            return self.dm.standard_df, self.dm.clinical_df
        
        else:
            if additional_population_info:
                return self.dm.get_additional_pop_info_df('raw'), self.dm.clinical_df
            return self.dm.raw_df, self.dm.clinical_df
        
    def get_json(self):
        variables = (self.dataset_id, self.transcript)
        return self.request_gnomad(variables)
    



class VariantSearch(Search):

    # variant_id: chromosome-position-original_nucleotide-variant
    #    example: 4-1002747-G-A 
    def __init__(self, dataset_version:int, variant_id: str):
        from pynoma.Queries import variant_search, variant_search_variables
        super().__init__(dataset_version, variant_search, variant_search_variables)
        self.variant_id = variant_id


    def get_json(self):
        variables = (self.dataset_id, self.variant_id)
        return self.request_gnomad(variables)


    def get_data(self, raw=False):
        json_data = self.get_json()
        if not json_data['data']['variant']:
            Logger.variant_not_found(self.variant_id)
            return (None, None)

        if raw:
            return json_data, None
        else:
            self.dm = DataManager(json_data, self.dataset_id, variant_search=True)
            return self.dm.standard_df, self.dm.variant_metadata
