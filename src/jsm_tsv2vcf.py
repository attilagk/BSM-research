#! /usr/bin/env python3

import csv
import math


def jsm_edit_records(jsm_tsv, out_tsv = "/dev/stdout", max_nrows = math.inf):
    with open(jsm_tsv, "r") as in_file, open(out_tsv, "w") as out_file:
        csv_reader = csv.reader(in_file, delimiter = "\t")
        csv_writer = csv.writer(out_file, delimiter = "\t")
        counter = 0
        m = []
        format_id = "counts_a:counts_b"
        for r in csv_reader:
            if counter > max_nrows:
                break
            d = ["."] # dot
            rr = r[0:2] + d + r[2:4] + d * 2
            normal = [ ":".join(str(x) for x in r[4:6]) ]
            tumor = [ ":".join(str(x) for x in r[6:8]) ]
            info = r[8:17]
            if counter == 0:
                info_id = info
            else:
                info_tagged = ";".join([ "=".join(x) for x in zip(info_id, info)])
                row = rr + [info_tagged] + [format_id] + normal + tumor
                csv_writer.writerow(row)
                m.append(row)
            counter += 1
        return(m)

if __name__ == "__main__":
    jsm_edit_records("/dev/stdin", "/dev/stdout")
