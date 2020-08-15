import pandas as pd

class DataManager:

    def __init__(self, json_data):
        self.json_data = json_data
        self.raw_df = None
        self.clinical_df = None
        self.standard_df = None
        self._process_raw_json()


    def _process_raw_json(self):
        clinical_var = self.json_data['data']['region']['clinvar_variants']
        variants = self.json_data['data']['region']['variants']
        self.raw_df = pd.DataFrame(variants)
        self.clinical_df = pd.DataFrame(clinical_var)
        return

    
    def process_standard_dataframe(self):
        self._process_raw_df()
        return

    
    def process_additional_pop_info(self):
        return


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
        return
        

    def _explicit_allele_informations(self, df):

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
                n_homs += population['ac_hom']
            
            num_homozygotes.append(n_homs)
            
        df['Allele Count'] = allele_count
        df['Allele Number'] = allele_number
        df['Allele Frequency'] = allele_freq
        df['Number of Homozygotes'] = num_homozygotes
        return df


