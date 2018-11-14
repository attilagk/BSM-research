#! /usr/bin/env python3

import csv
import math


def jsm_edit_records(jsm_tsv, out_tsv = "/dev/stdout", max_nrows = 5):
    with open(jsm_tsv, "r") as in_file, open(out_tsv, "w") as out_file:
        csv_reader = csv.reader(in_file, delimiter = "\t")
        csv_writer = csv.writer(out_file, delimiter = "\t")
        counter = 0
        for r in csv_reader:
            counter += 1
            if counter > max_nrows:
                break
            d = ["."] # dot
            rr = r[0:2] + d + r[2:4] + d * 2
            normal = [ ":".join(str(x) for x in r[4:6]) ]
            tumor = [ ":".join(str(x) for x in r[6:8]) ]
            csv_writer.writerow(rr + normal + tumor)
