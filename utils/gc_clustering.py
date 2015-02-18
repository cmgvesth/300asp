import csv
import numpy as np
genomes= dict()
with open('gc_processed.csv') as gc_data:
    for line in csv.reader(gc_data):
        genomes.setdefault(line[0],[]).append(line[1])

new = pd.DataFrame.from_dict(genome, orient='index', dtype=float)

new.to_csv('new')