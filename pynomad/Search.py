from requests import post

class Search:

    def __init__(self):
        self.end_point = "https://gnomad.broadinstitute.org/api/"

    
    def request_gnomad(self, query, variables):
        response = post(self.end_point, data={'query': query, 'variables': variables}, timeout=None)
        return response.json()


class RegionSearch(Search):

    def __init__(self, chromosome, start_position, end_position):
        Search.__init__(self)
        self.chromosome = str(chromosome)
        self.start = str(start_position)
        self.end = str(end_position)


class VariantSearch(Search):

    # variant_id: chromosome-position-original_nucleotide-variant
    #    example: 4-1002747-G-A 
    def __init__(self, variant_id: str):
        Search.__init__(self)
        self.variant_id = variant_id


class GeneSearch(Search):

    def __init__(self, gene: str):
        Search.__init__(self)
        self.gene = gene