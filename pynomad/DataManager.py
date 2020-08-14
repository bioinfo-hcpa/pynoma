import pandas as pd

class DataManager:

    def __init__(self):
        self.something = 1


    @classmethod
    def process_raw_json(self, received_json):
        clinical_var = received_json['data']['region']['clinvar_variants']
        variants = received_json['data']['region']['variants']
        return pd.DataFrame(variants), pd.DataFrame(clinical_var)


    @classmethod
    def get_raw_dataframes(self, received_json):
        return DataManager.process_raw_json(received_json)


    @classmethod
    def get_standard_dataframe(self, received_json):
        raw_df, clinical_df = DataManager.process_raw_json(received_json)
        standard_df = DataManager.process_raw_df(raw_df)
        return standard_df, clinical_df


    @classmethod
    def process_raw_df(self, raw_df):

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

        df_renamed = raw_df.rename(columns=renamed_cols)
        df_final = DataManager.explicit_allele_informations(df_renamed)
        return df_final.loc[:, standard_cols]
        

    @classmethod
    def explicit_allele_informations(self, df):

        allele_count = []
        allele_number = []
        allele_freq = []
        num_homozygotes = []

        for variant in df['genome']:
            allele_count.append(variant['ac'])
            allele_number.append(variant['an'])
            allele_freq.append("{:e}".format(variant['af']))
            
            n_homs = 0
            for population in variant['populations']:
                n_hom += population['ac_hom']
            
            num_homozygotes.append(n_homs)
            
        df['Allele Count'] = allele_count
        df['Allele Number'] = allele_number
        df['Allele Frequency'] = allele_freq
        df['Number of Homozygotes'] = num_homozygotes
        return df


