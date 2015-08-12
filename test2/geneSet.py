import argparse

def main():
    parser= argparse.ArgumentParser(description="Get file")
    parser.add_argument("file")
    args= parser.parse_args()
    f= open(args.file, "r")
    data= f.readlines()
    f.close()
    genes=[]
    for line in data:
        tokens = line.split(".")
        genes.append(tokens[0])

    outfile = open(args.file+"_unique.txt", "w")
    geneSet = set(genes)
    for g in geneSet:
        outfile.write(g+"\n")
    outfile.close()

main()
