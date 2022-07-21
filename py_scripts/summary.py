import os
import sys
from pathlib import Path
import pandas as pd

walk_dir = sys.argv[1]
i = 0
out = open("miRg3.0_summary.txt")
out.write("Sample name(s)\tTotal Input Reads\tTrimmed Reads (all)\tTrimmed Reads (unique)\tAll miRNA Reads\tFiltered miRNA Reads\tUnique miRNAs\tHairpin miRNAs\tmature tRNA Reads\tprimary tRNA Reads\tsnoRNA Reads\trRNA Reads\tncRNA others\tmRNA Reads\tRemaining Reads\n")
for root, subdirs, files in os.walk(walk_dir):
    for subdir in subdirs:
        if not "tRFs" in subdir:
            curDir = Path(Path(root).absolute()/subdir)
            miRC = str(Path(curDir/"annotation.report.csv"))
            with open(miRC, "r") as inf:
                rls = inf.readlines()
                for i in rls:
                    i = i.strip()
                    if "Sample name" not in i:
                        j = i.split(",")
                        i = "\t".join(j)
                        out.write(i + "\n")
