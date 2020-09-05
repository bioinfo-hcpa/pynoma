from requests import post
from pynomad.DataManager import DataManager

class Search:

    # dataset_version: either 3 or 2
    def __init__(self, dataset_version:int, query, query_variables):
        self.end_point = "https://gnomad.broadinstitute.org/api/"
        self.query = query
        self.query_vars = query_variables
        self.dataset_id = self.get_dataset_id(dataset_version)

        self.dm = None   # attribute holding DataManager object

    
    # variables: a tuple (var1, var2, ..., varn)
    def request_gnomad(self, variables):
        variables = self.query_vars % variables
        response = post(self.end_point, data={'query': self.query, 'variables': variables}, timeout=None)
        return response.json()

    
    def get_dataset_id(self, version):
        if version == 3:
            return 'gnomad_r3'
        elif version == 2:
            return 'gnomad_r2_1'
        else:
            raise Exception('Invalid dataset version. Choose either 2 or 3.')
        return


class RegionSearch(Search):

    # dataset_version: either 3 or 2
    def __init__(self, dataset_version:int, chromosome, start_position, end_position):

        from pynomad.Queries import in_region, in_region_variables
        super().__init__(dataset_version, in_region, in_region_variables)

        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)


    def get_json(self):
        variables = (self.chromosome, self.dataset_id, self.start, self.end)
        return self.request_gnomad(variables)


    # If standard is set to False, it will return everything without processing 
    # farther than {json to pandas DF}
    # If additional_population_info is set to True, 9 additional columns will 
    # be added to the returned dataframe (the 9 population allele frequency for each variant)
    def get_data(self, standard=True, additional_population_info=False):

        json_data = self.get_json()
        self.dm = DataManager(json_data)
        
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

        from pynomad.Queries import gene_id, gene_id_variables
        super().__init__(dataset_version, gene_id, gene_id_variables)
        
        self.gene = gene
        self.gene_ens_id = None
        self.get_ensembl_id(gene_id_variables)

        from pynomad.Queries import gene_search, gene_search_variables
        self.query = gene_search
        self.query_vars = gene_search_variables
        
        self.chromosome = None
        self.start = None
        self.end = None
        self.get_gene_information()

        self.region_search = RegionSearch(dataset_version, self.chromosome, 
                                            self.start, self.end)


    def get_ensembl_id(self, query_variables):
        json_data = self.request_gnomad((self.dataset_id, self.gene))
        str_ensg = json_data['data']['searchResults'][0]['value']
        self.gene_ens_id = str_ensg.split('/')[-1].split('?')[0]
        return


    def get_gene_information(self):
        gene_info = self.request_gnomad(self.gene_ens_id)
        self.chromosome = gene_info['data']['gene']['chrom']
        self.start = gene_info['data']['gene']['start']
        self.end = gene_info['data']['gene']['stop']
        return


    def get_data(self, standard=True, additional_population_info=False):
        dataframes = self.region_search.get_data(standard, additional_population_info)
        self.dm = self.region_search.dm
        return dataframes



class VariantSearch(Search):

    # variant_id: chromosome-position-original_nucleotide-variant
    #    example: 4-1002747-G-A 
    def __init__(self, dataset_version:int, variant_id: str):
        from pynomad.Queries import variant_search, variant_search_variables
        super().__init__(dataset_version, variant_search, variant_search_variables)
        self.variant_id = variant_id


    def get_json(self):
        variables = (self.dataset_id, self.variant_id)
        return self.request_gnomad(variables)