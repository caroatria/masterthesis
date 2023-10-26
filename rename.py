#!/usr/bin/env python3

import sys

def convToDict(info):
    init = iter(info.split(" "))
    dic = dict(zip(init, init))
    return dic

def parseRec(line):
    fields = line.strip().replace(";", "").split("\t")
    info = convToDict(fields[8])
    return(fields[0:8], info)

if __name__ == "__main__":
    fileName = sys.argv[1]
    names = dict()
    counter = 1;
    with open(fileName, "r") as fh:
        for line in fh:
            if "#" in line:
                continue
            f, i = parseRec(line)
            gid = ""
            tid = ""
            if i["gene_id"] not in names:
                gid = "SUB2.g" + str(counter)
                counter += 1
                names[i["gene_id"]] = { "gene_id" : gid }
            if i["transcript_id"] not in names[i["gene_id"]]:
                gid = names[i["gene_id"]]["gene_id"]
                tid = gid + ".t" + str(len(names[i["gene_id"]]))
                names[i["gene_id"]][i["transcript_id"]] = tid 
            gid = names[i["gene_id"]]["gene_id"]
            tid = names[i["gene_id"]][i["transcript_id"]]
            outLine = "\t".join(f) + "\tgene_id \"" + gid + "\"; transcript_id \"" + tid + "\";"
            print(outLine)
