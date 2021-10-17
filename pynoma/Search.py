from requests import post
from pynoma.DataManager import DataManager

class Search:

    # dataset_version: either 3 or 2
    def __init__(self, dataset_version:int, query, query_variables):
        self.end_point = "https://gnomad.broadinstitute.org/api/"
        self.query = query
        self.query_vars = query_variables
        self.reference_genome = None
        self.dataset_id = self.get_dataset_id(dataset_version)
        

        self.dm = None   # attribute holding DataManager object

    
    # variables: a tuple (var1, var2, ..., varn)
    def request_gnomad(self, variables):
        variables = self.query_vars % variables
        response = post(self.end_point, data={'query': self.query, 'variables': variables}, timeout=None)
        return response.json()

    
    def get_dataset_id(self, version):
        if version == 3:
            self.reference_genome = "GRCh38"
            return "gnomad_r3"
        elif version == 2:
            self.reference_genome = "GRCh37"
            return "gnomad_r2_1"
        else:
            raise Exception('Invalid dataset version. Choose either 2 or 3.')
        return


class RegionSearch(Search):

    # dataset_version: either 3 or 2
    def __init__(self, dataset_version:int, chromosome, start_position, end_position):

        from pynoma.Queries import in_region_v3, in_region_v2, in_region_variables
        if dataset_version == 2:
            in_region = in_region_v2
        else:
            in_region = in_region_v3
        super().__init__(dataset_version, in_region, in_region_variables)

        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)


    def get_json(self):
        variables = (self.chromosome, self.dataset_id, self.reference_genome, self.start, self.end)
        return self.request_gnomad(variables)


    # If standard is set to False, it will return everything without processing 
    # farther than {json to pandas DF}
    # If additional_population_info is set to True, 9 additional columns will 
    # be added to the returned dataframe (the 9 population allele frequency for each variant)
    def get_data(self, standard=True, additional_population_info=False):

        json_data = self.get_json()
        if not json_data['data']['region']['variants']:
            print("No variants found.")
            return (None, None)

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




class GeneSearch(Search):

    def __init__(self, dataset_version:int, gene: str):

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
            print("No gene found with given name.")
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
            print("No variants found.")
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
            print("No variants found for given transcript.")
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
            print("Variant not found.")
            return (None, None)

        if raw:
            return json_data, None
        else:
            self.dm = DataManager(json_data, self.dataset_id, variant_search=True)
            return self.dm.standard_df, self.dm.variant_metadata
