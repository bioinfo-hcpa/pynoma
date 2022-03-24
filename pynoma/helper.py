from pynoma.Logger import Logger
from matplotlib import style
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import pandas as pd

from random import uniform
from time import sleep


# path: "/my/saving/path/plot.png"
def annotation_barplot(df, path=None):
    style.use('fivethirtyeight')
    sns.set_palette("hls", len(df['Annotation'].value_counts().index))
    ax = sns.countplot(y="Annotation", data=df, order=df['Annotation'].value_counts().index)
    plt.tight_layout()
    
    if path:
        plt.savefig(path)

    plt.show()
    return ax


# search_objects: a list of Search objects different from VariantSearch, i.e,
# a list of GeneSearch, RegionSearch and/or TranscriptSearch objects
def batch_search(search_objects, standard=True, additional_population_info=False, verbose=True):
    datasets=[]
    total_searches = len(search_objects)
    for i, obj in enumerate(search_objects):
        if verbose:
            Logger.batch_searching(i+1, total_searches)
        sleep(uniform(1,5))
        try:
            obj_df, _ = obj.get_data(standard=standard, additional_population_info=additional_population_info) 
            if isinstance(obj_df, pd.DataFrame):
                datasets.append(obj_df)
        except Exception as e:
            if type(e).__name__ == 'KeyError':
                sleep(30)
                obj_df, _ = obj.get_data(standard=standard, additional_population_info=additional_population_info) 
                if isinstance(obj_df, pd.DataFrame):
                    datasets.append(obj_df)
            else:
                raise(e)
                
    if len(datasets) == 0:
        return None

    return pd.concat(datasets)#.fillna(0)