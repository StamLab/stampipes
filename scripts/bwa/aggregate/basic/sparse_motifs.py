import sys
import numpy as np
from sklearn.datasets import dump_svmlight_file

clusternames = sys.argv[1]
hotspots = sys.argv[2]

# hardcode output names
outfile = "hs_motifs_svmlight.txt"
outfile_rows = "hs_motifs_svmlight.rows.txt"
outfile_cols = "hs_motifs_svmlight.cols.txt"

with open(clusternames) as file:
    fimos = [line.strip() for line in file]

rows = []
hotspot_names = []
with open(hotspots) as file:
    for line in file:
        line_split = line.strip().split()
        hotspot_names.append(line_split[0] + "_" + line_split[1] + "_" + line_split[2])
        if len(line_split) == 10:
            motifs = line_split[9].split(";")
            newrow = [0] * len(fimos)
            for motif in motifs:
                index = [n for n, l in enumerate(fimos) if l == motif]
                newrow[index[0]] = 1
        else:
            newrow = [0] * len(fimos)
        rows.append(newrow)
labels = range(0, len(hotspot_names))

# write
dump_svmlight_file(rows, y=labels, f=outfile, zero_based=True)
file1 = open(outfile_rows, "w")
for item in hotspot_names:
    file1.write("%s\n" % item)
file2 = open(outfile_cols, "w")
for item in fimos:
    file2.write("%s\n" % item)
