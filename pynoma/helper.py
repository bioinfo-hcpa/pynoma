from matplotlib import style
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


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