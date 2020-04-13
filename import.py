# Created by Ivan Munoz-Gutierrez
# April 12, 2020
# Program to import csv features table into a db file for using SQL

import sqlite3
import csv
import sys

# checking correct usage of the Program
if len(sys.argv) != 2:
    sys.exit("usage: python import.py table.csv")

# creating a databse file and a connection to it
conn = sqlite3.connect("taxon.db")

# creating a cursor to execute sql
c = conn.cursor()

# create a features table in the database file
c.execute("""CREATE TABLE features (
             accession TEXT,
             organism TEXT,
             taxon INTEGER,
             seq_num INTEGER
             )""")

# we have to commit to save table
conn.commit()

# opening the csv features table
with open(sys.argv[1], "r") as infile:
    # using DictReader to facilitate the work
    reader = csv.DictReader(infile)

    # creating a features dictionary in memory
    features = {}
    counter = 1
    for row in reader:
        features[counter - 1] = {"accession": row["accession"],
                                 "organism": row["organism"],
                                 "taxon": int(row["taxon"]),
                                 "seq_num": counter
                                 }
        counter += 1
print("csv features table saved in memory as dictionary")

# saving features dictionary in database
for i in range(len(features)):
    c.execute("INSERT INTO features VALUES(?, ?, ?, ?)",
              (features[i]["accession"],
               features[i]["organism"],
               features[i]["taxon"],
               features[i]["seq_num"])
              )
    conn.commit()

conn.close()

print("csv features table saved in a db file")

sys.exit(0)
