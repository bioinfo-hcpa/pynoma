from requests import post
from pynomad.DataManager import DataManager

class Search:

    def __init__(self, dataset_id, query, query_variables):
        self.end_point = "https://gnomad.broadinstitute.org/api/"
        self.query = query
        self.query_vars = query_variables

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
            return 'gnomad_r2'
        else:
            raise('Invalid dataset version. Choose either 2 or 3.')
        return


class RegionSearch(Search):

    # dataset_version: either 3 or 2
    def __init__(self, dataset_version:int, chromosome, start_position, end_position):

        from pynomad.Queries import in_region, in_region_variables
        super().__init__(self, in_region, in_region_variables)

        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)
        
        self.dataset_id = super().get_dataset_id(dataset_version)


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



class VariantSearch(Search):

    # variant_id: chromosome-position-original_nucleotide-variant
    #    example: 4-1002747-G-A 
    def __init__(self, variant_id: str):
        super().__init__(self, "", "")
        self.variant_id = variant_id


class GeneSearch(Search):

    def __init__(self, gene: str):
        super().__init__(self, "", "")
        self.gene = gene