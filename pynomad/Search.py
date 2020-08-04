from requests import post

class Search:

    def __init__(self):
        self.end_point = "https://gnomad.broadinstitute.org/api/"

    
    def request_gnomad(self, query, variables):
        response = post(self.end_point, data={'query': query, 'variables': variables}, timeout=None)
        return response.json()


class RegionSearch(Search):

    def __init__(self):
        Search.__init__(self)


class VariantSearch(Search):

    def __init__(self):
        Search.__init__(self)


class GeneSearch(Search):

    def __init__(self):
        Search.__init__(self)