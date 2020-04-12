# Created by Ivan Munoz-Gutierrez
# April 11, 2020
# Program to fetch information from a collection of sequences from GeneBank.
# This program uses the history feature and the WebEnv session cookie to
# download large data in batches.

# importing Biopython modules
from Bio import Entrez
from Bio import SeqIO
import csv
import sys

# provide email address to GeneBank
Entrez.email = "email@address.com"

# using esearch to find the information, in this case we need
# data from the sequences in the BioProject 43021
# also we need to implement usehistory to retrieve large amount of sequences
search_handle = Entrez.esearch(db="nucleotide", term="43021[BioProject]",
                               usehistory="y")

# copying the information in computer memory
search_results = Entrez.read(search_handle)
search_handle.close()

# counting the number of results (sequences)
count = int(search_results["Count"])
print(count)

# copying cookie and query from history to keep track of our batch fetching
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

# number of sequences to be requested by batch.
# a batch of 500 is the max that we can request
batch_size = 500

# number to keep track of sequences, it is important in case the conextion
# to NCBI is interrupted so we can know where to continue downloading
seq_counter = 1

# opening our results file to write fetched data in csv format
with open("results.csv", "w") as results:
    writer = csv.writer(results)

    # field names or headers in the csv table
    fields = ["accession", "organism", "taxon", "seq_num"]

    # writing headers
    writer.writerow(fields)

    # fetching the information from Genebank by batches
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("going to download record %i to %i" % (start + 1, end))

        # getting the batch information
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc"
        )

        # parsing throw the information fetched
        for seq_record in SeqIO.parse(fetch_handle, "gb"):
            # dictionary features of a sequence is saved in feature
            feature = seq_record.features
            # looping throw dictionary feature
            for i in feature:
                # extracting sequnence id, i.e. accession number
                field0 = seq_record.id
                # looking for section source in feature
                if i.type == "source":
                    # creating a dictionary of the qualifiers from source
                    dictionary = dict(i.qualifiers)
                    # looking for organism name
                    if "organism" in dictionary:
                        organism = dictionary.get("organism")
                        field1 = organism[0]
                    else:
                        field1 = None
                    # looking for taxon number
                    if "db_xref" in dictionary:
                        taxones = dictionary.get("db_xref")
                        for k in range(len(taxones)):
                            if ":" in taxones[k]:
                                taxones_list = taxones[k].split(":")
                                if "taxon" in taxones_list:
                                    field2 = taxones_list[1]
                                    break
                    else:
                        field2 = None
                        break
                else:
                    break
            # keeping track of the number of sequence saved
            field3 = seq_counter
            seq_counter += 1
            fields = [field0, field1, field2, field3]
            # saving the retrived data in the csv file
            writer.writerow(fields)
        fetch_handle.close()

# if everything was OK and done print Done and exit the program
print("Done")
sys.exit(0)
