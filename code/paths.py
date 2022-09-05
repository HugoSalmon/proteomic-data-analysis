




from pathlib import Path
import os





source = Path().resolve().parent


for dir_ in ["data"]:
    if not os.path.exists(Path(source, dir_)):
        raise Exception("No " + dir_ + " directory in the proteomic_project directory")

for dir_ in ["figures", "save", "results"]:
    if not os.path.exists(Path(source, dir_)):
        print("No " + dir_ + " directory in the proteomic-data-analysis directory." \
                        " Creating one.")
    
        os.mkdir(Path(source, dir_))


datapath = Path(source, "data")
savepath = Path(source, "save")
figpath = Path(source, "figures") 
codepath = Path(source, "code")
resultspath = Path(source, "results")