# Created by Ivan Munoz-Gutierrez
# April 13, 2020
# Program to append several gb files with features
# The infile(s) must be in cvs format
# usage: python append_features_batches.py infile0.csv infile1.cvs infile(n-1).csv

import sys
import csv


# funtion to open input csv infile(s) and save it in memory
def saving(input_argument):
    with open(input_argument, "r") as csvinfile:
        reader = csv.DictReader(csvinfile)

        # copying data into computer memory
        features = {}
        counter = 0
        for row in reader:
            features[counter] = {"accession": row.get("accession"),
                                 "organism": row.get("organism"),
                                 "taxon": int(row.get("taxon")),
                                 }
            counter += 1
    return features


# saving the data in outfile
# looping throw argv to get the csvinfile(s)
counter = 1
for argument in range(len(sys.argv) - 1):
    # save cvsinfile in memory
    features = saving(sys.argv[argument + 1])

    # saving the data of the csvinfile
    with open("appended_results.csv", "a") as csvoutfile:
        fieldnames = ["accession", "organism", "taxon", "seq_num"]
        writer = csv.DictWriter(csvoutfile, fieldnames=fieldnames)

        # writing header
        if argument == 0:
            writer.writeheader()

        # writing cvsfile
        for i in range(len(features)):
            writer.writerow({"accession": features[i]["accession"],
                             "organism": features[i]["organism"],
                             "taxon": features[i]["taxon"],
                             "seq_num": counter
                             })
            counter += 1

print("done")
sys.exit(0)
