# created by Ivan Munoz-Gutierrez
# April 16, 2020
# program to remove rows from a pfam csv database
# it compares the organism in the row with a list of organims from a
# bioproject, if the organism is present in the bioproject the row is left
# additionally, repeated rows that have the same organim are removed

import csv
import sys

# checking the correct usage of the program
if len(sys.argv) != 3:
    sys.exit("Usage: python bioproject_file.csv pfam_file.csv")

# opening and saving the organims of the bioproject csv file
# in a list of dictionaries
bioproject = []
with open(sys.argv[1], "r") as infile_bioproject:
    reader = csv.DictReader(infile_bioproject)

    # copying the infile into a dictionary
    for row in reader:
        bacteria = {"organism": row["organism"]}
        bioproject.append(bacteria)

# opening pfam file and creating a list of dictionaries
list_pfam = []
with open(sys.argv[2], "r") as infile_pfam:
    reader = csv.DictReader(infile_pfam)

    # saving the csv table in list_pfam
    for row in reader:
        list_pfam.append(row)

# sorting the list_pfam by "Organims"
def myFunc(e):
    return e["Organism"]

list_pfam.sort(key=myFunc)

# deleting duplicates in list_pfam and making a dictionary
dict_pfam = {}
counter = 0
for i in range(len(list_pfam)):
    # saving the first row
    if counter == 0:
        dict_pfam[counter] = list_pfam[i]
        counter += 1
    else:
        strain_A = list_pfam[i - 1]["Organism"].split(" ")
        # adding a second index in case the name of Organims is not clear.
        # sometimes the Organism filed says only bacterium
        if len(strain_A) == 1:
            strain_A.append(None)
        strain_B = list_pfam[i]["Organism"].split(" ")
        if len(strain_B) == 1:
            strain_B.append(None)
        # selecting row with bacteria from the same genus but different species
        if strain_A[0] == strain_B[0] and strain_A[1] != strain_B[1]:
            dict_pfam[counter] = list_pfam[i]
            counter += 1
        # selecting bacteria from different genus
        elif strain_A[0] != strain_B[0]:
            dict_pfam[counter] = list_pfam[i]
            counter += 1

# getting the headers of the pfam file
fieldnames = []
with open(sys.argv[2], "r") as infile_pfam:
    reader = csv.reader(infile_pfam)

    # getting the headers of the pfram file
    fieldnames = next(reader)

# saving the relevant information from the pfam file into a new csv outfile
with open("clean_pfam.csv", "w") as outfile_pfam:
    writer = csv.DictWriter(outfile_pfam, fieldnames=fieldnames)

    # writing header
    writer.writeheader()

    for row in range(len(dict_pfam)):
        # copying the genus and species of row
        bacterium_pfam = dict_pfam[row]["Organism"].split(" ")
        # checking if bacterium_pfam is in the bioproject
        for i in range(len(bioproject)):
            bacterium_bioproject = bioproject[i]["organism"].split(" ")
            if bacterium_pfam[0] == bacterium_bioproject[0] and bacterium_pfam[1] == bacterium_bioproject[1]:
                writer.writerow(dict_pfam[row])

print("Done")
sys.exit(0)
