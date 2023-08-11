import importlib_resources
import pandas as pd

def load_markers():

    my_resources = importlib_resources.files("sctools.tools") / "src"
    return(pd.read_csv((my_resources / "gene_markers.csv")))
