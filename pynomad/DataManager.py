import pandas as pd
from copy import deepcopy


POPULATION_ID_MAP = {
            'AFR': 'African',
            'AMI': 'Amish',
            'AMR': 'Latino',
            'ASJ': 'Ashkenazi Jewish',
            'EAS': 'East Asian',
            'FIN': 'European (Finnish)',
            'NFE': 'European (non-Finnish)',
            'OTH': 'Other',
            'SAS': 'South Asian'
        }



class DataManager:

    def __init__(self, json_data, variant_search=False):
        self.json_data = json_data
        self.variant_search = variant_search

        self.raw_df = None
        self.clinical_df = None
        self.standard_df = None
        self.variant_metadata = None  # Used in variant searches
        self._process_raw_json()


    def _process_raw_json(self):
        if self.variant_search:
            self.__process_variant_search_data()
        
        else:
            clinical_var = self.json_data['data']['region']['clinvar_variants']
            variants = self.json_data['data']['region']['variants']
            self.raw_df = pd.DataFrame(variants)
            self.clinical_df = pd.DataFrame(clinical_var)
        return

    
    def process_standard_dataframe(self):
        self._process_raw_df()
        return

    
    # dataframe_type: either 'standard' or 'raw' (it is not binary because 
    # maybe more output dataframe types will be added in the future)
    def get_additional_pop_info_df(self, dataframe_type:str):

        populations_freq_column = self._process_populations_frequency()
        if dataframe_type == 'standard':
            return self._add_populations_freq_column(populations_freq_column, self.standard_df)
        elif dataframe_type == 'raw':
            return self._add_populations_freq_column(populations_freq_column, self.raw_df)
        else:
            raise("Dataframe type should be either 'raw' or 'standard'")


    
    def _process_populations_frequency(self):
        
        populations_freq_column = {}
        for pop_id in POPULATION_ID_MAP:
            populations_freq_column[pop_id] = []
    
        for variant in self.raw_df['genome']:
            for population in variant['populations']:
                pop_id = population['id']
                pop_allele_freq = population['ac'] / population['an']
                populations_freq_column[pop_id].append("{:e}".format(pop_allele_freq))

        return populations_freq_column


    def _add_populations_freq_column(self, populations_freq_column, df):
        df_copy = deepcopy(df)
        for pop_id, pop_freqs in populations_freq_column.items():
            df_copy[POPULATION_ID_MAP[pop_id]] = pop_freqs
        return df_copy


    def _process_raw_df(self):

        renamed_cols = {
                            'variant_id': 'Variant ID',
                            'gene_symbol': 'Gene',
                            'hgvs': 'Consequence',
                            'consequence': 'Annotation',
                            'flags': 'Flags'
                        }
        
        standard_cols = [
                    'Variant ID', 'Gene', 'Consequence', 
                    'Annotation', 'Flags', 'Allele Count',
                    'Allele Number', 'Allele Frequency',
                    'Number of Homozygotes'
                ]

        df_renamed = self.raw_df.rename(columns=renamed_cols)
        df_final = self._explicit_allele_informations(df_renamed)
        self.standard_df = df_final.loc[:, standard_cols]
        self._add_variant_columns()
        return
        

    def _explicit_allele_informations(self, df):

        allele_count = []
        allele_number = []
        allele_freq = []
        num_homozygotes = []

        for variant in df['genome']:
            allele_count.append(variant['ac'])
            allele_number.append(variant['an'])
            #allele_freq.append("{:e}".format(variant['af']))
            allele_freq.append(variant['af'])
            
            n_homs = 0
            for population in variant['populations']:
                n_homs += population['ac_hom']
            
            num_homozygotes.append(n_homs)
            
        df['Allele Count'] = allele_count
        df['Allele Number'] = allele_number
        df['Allele Frequency'] = allele_freq
        df['Number of Homozygotes'] = num_homozygotes
        return df


    def _add_variant_columns(self):

        chromosome_col = []
        location_col = []
        gen_reference_col = []
        gen_alternative_col = []

        for variant in self.standard_df['Variant ID']:
            info_pieces = variant.split('-')
            chromosome_col.append(info_pieces[0])
            location_col.append(info_pieces[1])
            gen_reference_col.append(info_pieces[2])
            gen_alternative_col.append(info_pieces[3])

        self.standard_df['Chromosome'] = chromosome_col
        self.standard_df['Location'] = location_col
        self.standard_df['Reference'] = gen_reference_col
        self.standard_df['Alternative'] = gen_alternative_col
        return
    

    def __process_variant_search_data(self):
        self.__build_variant_search_standard_df()
        self.__extract_variant_metadata()
        return
    

    def __build_variant_search_standard_df(self):
        
        reqdf = pd.json_normalize(self.json_data['data']['variant']['genome']['populations']).set_index('id')
        POPULATION_ID_MAP['FEMALE'] = 'Female'
        POPULATION_ID_MAP['MALE'] = 'Male'

        new_index = {}
        frequencies = []
        for row in reqdf.iterrows():
            row_id = row[0]
            new_row_name = ""
            for piece in row_id.split('_'):
                new_row_name += POPULATION_ID_MAP[piece] + " "
            new_row_name = new_row_name[0:-1]
            new_index[row_id] = new_row_name
            frequencies.append(row[1]['ac']/row[1]['an'])

        new_index['FEMALE'] = 'Total Female'
        new_index['MALE'] = 'Total Male'
        new_columns = {'ac': 'Allele Count', 'an': 'Allele Number',
                    'ac_hemi': 'Number of Hemizygotes', 'ac_hom': 'Number of Homozygotes'}

        df = reqdf.rename(columns=new_columns, index=new_index)

        total_row = df.loc['Total Female'] + df.loc['Total Male']
        total_row.name = 'Total'
        df = df.append([total_row])

        frequencies.append(df.loc['Total']['Allele Count'] / df.loc['Total']['Allele Number'])
        df['Allele Frequency'] = frequencies
        
        chromosome = self.json_data['data']['variant']['chrom']
        if chromosome != 'X' and chromosome != 'Y':
            del df['Number of Hemizygotes']

        self.standard_df = df
        return


    def __extract_variant_metadata(self):
        metadata = deepcopy(self.json_data)
        metadata['data']['variant']['genome'].pop('populations')
        self.variant_metadata = metadata['data']['variant']
        return
