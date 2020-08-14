from requests import post
from pynomad.DataManager import DataManager

class Search:

    def __init__(self, dataset_id, query, query_variables):
        self.end_point = "https://gnomad.broadinstitute.org/api/"
        self.query = query
        self.query_vars = query_variables

    
    # variables is a tuple (var1, var2, ..., varn)
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

        from Queries import in_region, in_region_variables
        super().__init__(self, in_region, in_region_variables)

        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)
        
        self.dataset_id = super().get_dataset_id(dataset_version)


    def get_json(self):
        variables = (self.chromosome, self.dataset_id, self.start, self.end)
        return self.request_gnomad(variables)


    def get_dataframe(self, standard=True, additional_population_info=False, clinical_dataframe=False):

        json_data = self.get_json()
        df = None
        clinical_df = None 

        if standard:
            df, clinical_df = DataManager.get_standard_dataframe(json_data)
        else:
            df, clinical_df = DataManager.get_raw_dataframes(json_data)
        
        if additional_population_info:
            df = df #TO-DO!!!

        if clinical_dataframe:
            return df, clinical_df
        else:
            return df


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