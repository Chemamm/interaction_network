import NetworkMaster
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',help='tsv file')
parser.add_argument('-o','--output',help='output label for filenames')
parser.add_argument('-n','--nclusters',help='number of clusters to generate using certain algorithms')
parser.add_argument('-s', '--source',help='Source node column name')
parser.add_argument('-t', '--target',help='Target node column name')
args = parser.parse_args()

G = NetworkMaster.network_import_directed(args.file,args.source,args.target)
NetworkMaster.network_stats(G, args.output + "_network_stats.tsv")
NetworkMaster.network_community_analysis(G, args.output)