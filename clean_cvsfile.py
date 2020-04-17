# created by Ivan Munoz-Gutierrez
# April 16, 2020
# program to remove rows that have the same organism

import csv
import sys

# checking the correct usage of the Program
if len(sys.argv) != 2:
    sys.exit("Usage: python clean_csvfile.csv")

# open csv infile and
# creating a list of dictionaries of features
list_features = []
with open(sys.argv[1], "r") as csvinfile:
    reader = csv.DictReader(csvinfile)

    # saving the csv table in list_features
    for row in reader:
        features = {"accession": row["accession"],
                    "organism": row["organism"],
                    "taxon": row["taxon"],
                    "seq_num": row["seq_num"]
                    }
        list_features.append(features)

# sorting the list_features by "organims"
def myFunc(e):
    return e["organism"]

list_features.sort(key=myFunc)

# deleting duplicated species from the list and
# making the dictionary features to export it as cvs outfile
features = {}
counter = 0
for i in range(len(list_features)):
    # saving features of the first organism
    if counter == 0:
        features[counter] = {"accession": list_features[i]["accession"],
                             "organism": list_features[i]["organism"],
                             "taxon": list_features[i]["taxon"],
                             "seq_num": list_features[i]["seq_num"]
                             }
        counter += 1
    # selecting bacteria from the same genus but different species
    else:
        strain_A = list_features[i - 1]["organism"].split(" ")
        strain_B = list_features[i]["organism"].split(" ")
        if strain_A[0] == strain_B[0] and strain_A[1] != strain_B[1]:
            features[counter] = {"accession": list_features[i]["accession"],
                                 "organism": list_features[i]["organism"],
                                 "taxon": list_features[i]["taxon"],
                                 "seq_num": list_features[i]["seq_num"]
                                 }
            counter += 1
        # selecting bacteria from different genus
        elif strain_A[0] != strain_B[0]:
            features[counter] = {"accession": list_features[i]["accession"],
                                 "organism": list_features[i]["organism"],
                                 "taxon": list_features[i]["taxon"],
                                 "seq_num": list_features[i]["seq_num"]
                                 }
            counter += 1

# copying the clean csv table into the csvoutfile
with open("clean_csvfile.csv", "w") as csvoutfile:
    fieldnames = ["accession", "organism", "taxon", "seq_num"]
    writer = csv.DictWriter(csvoutfile, fieldnames=fieldnames)

    # writing header
    writer.writeheader()

    # writing cvsfile
    counter = 0
    for i in range(len(features)):
        writer.writerow({"accession": features[i]["accession"],
                         "organism": features[i]["organism"],
                         "taxon": features[i]["taxon"],
                         "seq_num": counter
                         })
        counter += 1

print("done")
sys.exit(0)
