# Program to fetch information from Genebank
# Created by Ivan Munoz-Gutierrez
# April 5, 2020

from Bio import Entrez
from Bio import SeqIO
import csv
import sys

# checking the correct useage of the program
if len(sys.argv) != 2:
    sys.exit("usage: python fetch_aqc_MSP.py accession_list.txt")

# opening the infile.txt list
reader = open(sys.argv[1], 'r')

# creating a list of accession numbers in memory
accession_number = []
for row in reader:
    accession_number.append(row.replace('\n', ''))

# closing the infile.txt
reader.close()

# working with Genebank
# providing email address to Genebank
Entrez.email = "ivan.munoz.gutierrez@gmail.com"

# creating file to save the requested data from Genebank
with open("results.csv", "w") as results:

    writer = csv.writer(results)

    # field names
    fields = ["accession", "size", "topology", "organism", "strain",
              "isolation_source", "host", "plasmid", "country"]

    # writing headers
    writer.writerow(fields)

    # requesting the information from the Genebank
    for request in range(len(accession_number) - 1):
        field0 = accession_number[request + 1]
        print(field0)

        with Entrez.efetch(
            db="nucleotide", rettype="gb", retmode="text", id=accession_number[request + 1]
        ) as handle:
            record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

        # .seq is an object with the sequence itself
        field1 = len(record.seq)

        # .annotations is a dictionary of aditinonal information about the sequence as
        # topology, sequence_version, organims, references, etc.
        if "topology" in record.annotations:
            field2 = record.annotations["topology"]
        else:
            field2 = None

        # .features is a list of SeqFeatures objest with more structured information
        # about the features on a sequence
        feature = record.features

        for source in feature:
            # .type is only a description of the type of feature
            # that could be source, CDS, gene, etc.
            # in source we can find organism, strain, host, country, etc.
            if source.type == 'source':
                dictionary = dict(source.qualifiers)  # making a dictionary of qualifiers
                # .get gives a list
                if 'organism' in dictionary:
                    bacterium = dictionary.get('organism')
                    field3 = bacterium[0]
                else:
                    field3 = None

                if 'strain' in dictionary:
                    strain = dictionary.get('strain')
                    field4 = strain[0]
                else:
                    field4 = None

                if 'isolation_source' in dictionary:
                    isolation_source = dictionary.get('isolation_source')
                    field5 = isolation_source[0]
                else:
                    field5 = None

                if 'host' in dictionary:
                    host = dictionary.get('host')
                    field6 = host[0]
                else:
                    field6 = None

                if 'plasmid' in dictionary:
                    plasmid = dictionary.get('plasmid')
                    field7 = plasmid[0]
                else:
                    field7 = None

                if 'country' in dictionary:
                    country = dictionary.get('country')
                    field8 = country[0]
                else:
                    field8 = None

        fields = [field0, field1, field2, field3, field4, field5, field6, field7, field8]
        writer.writerow(fields)

print("Done")

sys.exit(0)
