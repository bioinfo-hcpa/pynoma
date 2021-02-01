from matplotlib import style
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import pandas as pd


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
def batch_search(search_objects, standard=True, additional_population_info=False):
    datasets=[]
    for obj in search_objects:
        obj_df, _ = obj.get_data(standard=standard, additional_population_info=additional_population_info) 
        if isinstance(obj_df, pd.DataFrame):
            datasets.append(obj_df)
    return pd.concat(datasets)#.fillna(0)