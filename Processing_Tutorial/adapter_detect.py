import subprocess
import os
import sys
from pathlib import Path
import pandas as pd

walk_dir = sys.argv[1]

for root, subdirs, files in os.walk(walk_dir):
	for ifin in files:
		diry = str(walk_dir)+"/"+str(ifin)
		command = "zgrep \"TGAGGTAGTAGGTTGTATAGTT\" " + str(diry) + " | head -5 "
		output = str(subprocess.check_output(command, shell=True))
		val = output.split("'")[1]
		val1 = val.split("\\n")
		for x in val1:
			if str(x) != '':   
				print(str(x) + "\t" + str(ifin))

