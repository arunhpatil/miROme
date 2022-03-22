import os
import sys
from pathlib import Path
import pandas as pd

walk_dir = sys.argv[1]

i = 0
for root, subdirs, files in os.walk(walk_dir):
    for subdir in subdirs:
        if not "tRFs" in subdir:
            i += 1
            curDir = Path(Path(root).absolute()/subdir)
            miRC = str(Path(curDir/"miR.RPM.csv"))
            if i == 1:
                d = pd.read_csv(miRC)
                d.set_index('miRNA',inplace = True)
            else:
                f = pd.read_csv(miRC)
                f.set_index('miRNA',inplace = True)
                d = d.join(f, how='outer',lsuffix='_l', rsuffix='_r')
d.to_csv("combined_RPM.csv")
