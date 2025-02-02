import NetworkMaster
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',help='tsv file')
parser.add_argument('-o','--output',help='output label for filenames')
parser.add_argument('-m','--map',help='Seed mapping file')
parser.add_argument('-e','--ematrix', default=False, help='Ematrix file')
parser.add_argument('-p','--poscon',help='Positional consensus file')
args = parser.parse_args()

fl = NetworkMaster.ClusteredNetwork(args.file)
fl.get_full_cluster(args.poscon, args.map,args.ematrix, args.output)
fl.get_enrichment(args.output, args.poscon)
fl.venn(args.poscon, args.output)