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
            'SAS': 'South Asian',
            'MID': 'Middle Eastern'
        }

SUBPOPULATION_ID_MAP = {
            'XX': 'XX',
            'XY': 'XY',
            'FEMALE': 'XX',
            'MALE': 'XY',

            'JPN': 'Japanese',
            'KOR': 'Korean',
            'OEA': 'Other',

            'BGR': 'Bulgarian',
            'EST': 'Estonian',
            'NWE': 'North-western',
            'SEU': 'Southern',
            'SWE': 'Swedish',
            'ONF': 'Other'
        }



class DataManager:

    def __init__(self, json_data, gnomad_version:str, variant_search=False, second_level_key='region'):
        self.json_data = json_data
        self.variant_search = variant_search
        self.gnomad_version = gnomad_version 
        self.second_level_key = second_level_key

        self.raw_df = None
        self.clinical_df = None
        self.standard_df = None
        self.variant_metadata = None  # Used in variant searches
        self._process_raw_json()


    def _process_raw_json(self):
        if self.variant_search:
            self.__process_variant_search_data()
        
        else:
            clinical_var = self.json_data['data'][self.second_level_key]['clinvar_variants']
            variants = self.json_data['data'][self.second_level_key]['variants']
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

        for row in self.raw_df.loc[:, ['genome','exome']].iterrows():
            
            genome_variant = row[1]['genome']
            exome_variant = row[1]['exome']

            if genome_variant:
                if exome_variant:   # thre's genome and exome data
                    for i in range(len(genome_variant['populations'])):
                        ac = exome_variant['populations'][i]['ac'] + genome_variant['populations'][i]['ac']
                        an = exome_variant['populations'][i]['an'] + genome_variant['populations'][i]['an']
                        pop_id = genome_variant['populations'][i]['id']
                        pop_allele_freq = ac/an if an else 0
                        populations_freq_column[pop_id.upper()].append("{:e}".format(pop_allele_freq))
                
                else:   # there's only genome data
                    for i in range(len(genome_variant['populations'])):
                        ac = genome_variant['populations'][i]['ac']
                        an = genome_variant['populations'][i]['an']
                        pop_id = genome_variant['populations'][i]['id']
                        pop_allele_freq = ac/an if an else 0
                        populations_freq_column[pop_id.upper()].append("{:e}".format(pop_allele_freq))
        

            else:  # there's only exome data
                for i in range(len(exome_variant['populations'])):
                    ac = exome_variant['populations'][i]['ac']
                    an = exome_variant['populations'][i]['an']
                    pop_id = exome_variant['populations'][i]['id']
                    pop_allele_freq = ac/an if an else 0
                    populations_freq_column[pop_id.upper()].append("{:e}".format(pop_allele_freq))

        return populations_freq_column


    def _add_populations_freq_column(self, populations_freq_column, df):
        df_copy = deepcopy(df)
        for pop_id, pop_freqs in populations_freq_column.items():
            if pop_freqs:   # If there is actually recorded frequencies for that population
                df_copy[POPULATION_ID_MAP[pop_id]] = pop_freqs
        return df_copy


    def _process_raw_df(self):

        renamed_cols = {
                            'variant_id': 'Variant ID',
                            'rsid': 'rsID',
                            'gene_symbol': 'Gene',
                            'hgvs': 'Consequence',
                            'consequence': 'Annotation',
                            'flags': 'Flags'
                        }
        
        standard_cols = [
                    'Variant ID', 'rsID', 'Gene', 'Consequence', 
                    'Annotation', 'Flags', 'Allele Count',
                    'Allele Number', 'Allele Frequency',
                    'Number of Homozygotes'
                ]

        

        df_renamed = self.raw_df.rename(columns=renamed_cols)
        df_final = self._explicit_allele_informations(df_renamed, standard_cols)
        self.standard_df = df_final.loc[:, standard_cols]
        self._add_variant_columns()
        return
        
    

    def _explicit_allele_informations(self, df, standard_cols):
        chromosome = df['Variant ID'][0][0]
        allele_count = []
        allele_number = []
        allele_freq = []
        num_homozygotes = []
        num_hemizygotes = []
        region = []

        for row in df.loc[:, ['genome','exome']].iterrows():
            
            genome_variant = row[1]['genome']
            exome_variant = row[1]['exome']
            n_homs = 0
            n_hemi = 0
            
            if genome_variant:
                n_homs, n_hemi = self._count_homos_hemis_variant_pops(n_homs, n_hemi, genome_variant)
                
                if exome_variant:   #there's genome and exome
                    n_homs, n_hemi = self._count_homos_hemis_variant_pops(n_homs, n_hemi, exome_variant)
                    variant_ac = genome_variant['ac']+exome_variant['ac']
                    variant_an = genome_variant['an']+exome_variant['an']
                    variant_af = genome_variant['af']+exome_variant['af']
                    region.append("Genome and Exome")
                    
                else:   #there's only genome
                    variant_ac = genome_variant['ac']
                    variant_an = genome_variant['an']
                    variant_af = genome_variant['af']
                    region.append("Genome")
                    
            else:   #there's only exome
                n_homs, n_hemi = self._count_homos_hemis_variant_pops(n_homs, n_hemi, exome_variant)
                variant_ac = exome_variant['ac']
                variant_an = exome_variant['an']
                variant_af = exome_variant['af']
                region.append("Exome")
        
            num_homozygotes.append(n_homs)
            num_hemizygotes.append(n_hemi)
            allele_count.append(variant_ac)
            allele_number.append(variant_an)
            allele_freq.append(variant_af)


        df['Allele Count'] = allele_count
        df['Allele Number'] = allele_number
        df['Allele Frequency'] = allele_freq
        df['Number of Homozygotes'] = num_homozygotes
        if (chromosome == 'X') or (chromosome == 'Y'):
            df['Number of Hemizygotes'] = num_hemizygotes
            standard_cols.append('Number of Hemizygotes')
        df['Source'] = region

        #if self.gnomad_version == 'gnomad_r2_1':
        standard_cols.append('Source')
        return df

    
    def _count_homos_hemis_variant_pops(self, n_homs, n_hemi, variant):
        for population in variant['populations']:
            n_homs += population['ac_hom']
            n_hemi += population['ac_hemi']
        return n_homs, n_hemi


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
        
        if self.json_data['data']['variant']['genome']:     
            if self.json_data['data']['variant']['exome']:  # genome and exome
                reqdf = pd.json_normalize(self.json_data['data']['variant']['genome']['populations']).set_index('id') + \
                        pd.json_normalize(self.json_data['data']['variant']['exome']['populations']).set_index('id')
            else:  # only genome
                reqdf = pd.json_normalize(self.json_data['data']['variant']['genome']['populations']).set_index('id')
        else:   # only exome
            reqdf = pd.json_normalize(self.json_data['data']['variant']['exome']['populations']).set_index('id')

        new_index = {}
        frequencies = []
        for row in reqdf.iterrows():
            row_id = row[0]
            new_row_name = ""
            
            for piece in row_id.split('_'):

                piece = piece.upper()
                try:
                    new_row_name += POPULATION_ID_MAP[piece] + " "
                except:
                    try:
                        new_row_name += SUBPOPULATION_ID_MAP[piece] + " "
                    except:
                        new_row_name += piece

            new_row_name = new_row_name[0:-1]
            new_index[row_id] = new_row_name
            
            frequency = 0
            if row[1]['an']:
                frequency = row[1]['ac']/row[1]['an']
            frequencies.append(frequency)
            

        new_index['XX'] = 'Total XX'
        new_index['XY'] = 'Total XY'
        new_index['FEMALE'] = 'Total XX'
        new_index['MALE'] = 'Total XY'
        new_columns = {'ac': 'Allele Count', 'an': 'Allele Number',
                    'ac_hemi': 'Number of Hemizygotes', 'ac_hom': 'Number of Homozygotes'}

        df = reqdf.rename(columns=new_columns, index=new_index)
        
        total_row = df.loc['Total XX'] + df.loc['Total XY']
        total_row.name = 'Total'
        df = df.append([total_row])

 
        frequencies.append(df.loc['Total']['Allele Count'] / df.loc['Total']['Allele Number'])
        df['Allele Frequency'] = frequencies
        df['Allele Frequency'] = df['Allele Frequency'].fillna(0)
        
        chromosome = self.json_data['data']['variant']['chrom']
        if chromosome != 'X' and chromosome != 'Y':
            del df['Number of Hemizygotes']

        self.standard_df = df
        return


    def __extract_variant_metadata(self):
        metadata = deepcopy(self.json_data)
        if metadata['data']['variant']['genome']:
            metadata['data']['variant']['genome'].pop('populations')
        if metadata['data']['variant']['exome']:
            metadata['data']['variant']['exome'].pop('populations')
        self.variant_metadata = metadata['data']['variant']
        return
