import importlib_resources
import pandas as pd

def load_markers():
    my_resources = importlib_resources.files("sctools.tools") / "src"
    #text = (my_resources / "gene_markers.csv").read_text()
    #template = importlib.resources.read_text("tools.src", "gene_markers.csv")
    return(pd.read_csv((my_resources / "gene_markers.csv")))
