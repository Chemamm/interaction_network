import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import community
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.community import k_clique_communities
from networkx.algorithms.community import label_propagation_communities
from networkx.algorithms.community import asyn_lpa_communities
from networkx.algorithms.community import greedy_modularity_communities
import stringdb
import venn
from scipy.stats import fisher_exact
from collections import Counter



class ClusteredNetwork:
    def __init__(self, file):
        with open(file, 'r') as fl:
            self.header = fl.readline()
            self.lines = fl.readlines()[1::]
        self.seed_dict, self.gene_dict = self.get_cluster()

    class Line:
        def __init__(self, line):
            self.fields = line.strip('\n').split(',')
            self.cluster = int(self.fields[0].strip('"'))
            self.node = self.fields[1].strip('"')

    def get_cluster(self):
        seed_dict = {}
        gene_dict = {}
        for line in self.lines:
            line = ClusteredNetwork.Line(line)
            if line.cluster not in seed_dict:
                seed_dict[line.cluster] = []
                gene_dict[line.cluster] = []
            if "seed" in line.node:
                seed_dict[line.cluster].append(line.node)
            else:
                gene = line.node
                if "." in gene:
                    gene = gene.split(".")[0]
                if "|" in gene:
                    gene = gene.split("|")[0]
                gene_dict[line.cluster].append(gene)
        return seed_dict, gene_dict

    def find_intersection_genes(self, posconsensus):
        intersection_cluster_dict = {}
        common_genes_dict = {}
        poscon = Positional_Consensus(posconsensus)

        for cluster, mirna_list in self.seed_dict.items():
            if mirna_list:
                first_mirna = mirna_list[0]
                cluster_genes_set = set(poscon.mirna_dict.get(first_mirna, []))
                #print(cluster_genes_set)

                # Iterate through the remaining miRNAs in the cluster
                for mirna in mirna_list[1:]:
                    # Find the intersection of genes with the current miRNA
                    cluster_genes_set.intersection_update(poscon.mirna_dict.get(mirna, []))
                    #print(cluster_genes_set)
                cluster_genes_set.intersection_update(self.gene_dict[cluster])
                # Store the result in the intersection_cluster_dict
                intersection_cluster_dict[cluster] = cluster_genes_set

                mirna_genes = [set(poscon.mirna_dict.get(mirna, [])) for mirna in mirna_list]

                # Find common genes in at least two lists
                common_genes = find_common_genes(mirna_genes, min_occurrence=2)
                common_genes_dict[cluster] = [x for x in common_genes if x in self.gene_dict[cluster]]


        return intersection_cluster_dict, common_genes_dict


    def get_full_cluster(self, posconsensus, mirna_map_file, ematrix, output):
        intersection_cluster_dict, common_genes = self.find_intersection_genes(posconsensus)
        if ematrix:
            ematrix = Ematrix(ematrix)
            ematrix.get_median()
            with open(output, 'wt') as out:
                out.write("Cluster\tSeeds\tNumSeeds\tmiRNAs\tMedian Abundance by MiRNA\tMean of Medians\tNumGenes"
                          "\tNumCommonGenes\tNumIntersectionGenes\tGenes\tCommonGenes\tIntersectionGenes\n")
                for cluster in self.seed_dict:
                    mirnas = seed_list_to_mirna(self.seed_dict[cluster], mirna_map_file)
                    subset_df = ematrix.df[ematrix.df['name'].isin(mirnas)]

                    # Calculate the median along the columns for each miRNA
                    median_values = subset_df.iloc[:, 1:].median(axis=1)
                    mean_of_median_values = median_values.mean()


                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cluster, ",".join(self.seed_dict[cluster]),
                                                                                len(self.seed_dict[cluster]),
                                                                            ",".join(mirnas),
                                                                                        ",".join([str(round(x,2)) for x in median_values]),
                                                                                        mean_of_median_values,
                                                                            len(self.gene_dict[cluster]),
                                                                            len(common_genes[cluster]),
                                                                            len(intersection_cluster_dict[cluster]),
                                                                ",".join(self.gene_dict[cluster]),
                                                                ",".join(common_genes[cluster]),
                                                                ",".join(intersection_cluster_dict[cluster])))
        else:
            with open(output, 'wt') as out:
                out.write("Cluster\tSeeds\tNumSeeds\tmiRNAs\tNumGenes"
                          "\tNumCommonGenes\tNumIntersectionGenes\tGenes\tCommonGenes\tIntersectionGenes\n")
                for cluster in self.seed_dict:
                    mirnas = seed_list_to_mirna(self.seed_dict[cluster], mirna_map_file)

                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cluster, ",".join(
                        self.seed_dict[cluster]),
                                                                                        len(self.seed_dict[cluster]),
                                                                                        ",".join(mirnas),
                                                                                        len(self.gene_dict[cluster]),
                                                                                        len(common_genes[cluster]),
                                                                                        len(intersection_cluster_dict[
                                                                                                cluster]),
                                                                                        ",".join(
                                                                                            self.gene_dict[cluster]),
                                                                                        ",".join(common_genes[cluster]),
                                                                                        ",".join(
                                                                                            intersection_cluster_dict[
                                                                                                cluster])))

    def get_enrichment(self, output, poscon=False):
        for cluster in self.gene_dict:
            if len(self.gene_dict[cluster]) >= 2:
                print('normal: cluster%i, len %i' %(cluster, len(self.gene_dict[cluster])))
                df = get_functional_enrichment(self.gene_dict[cluster], cluster)
                if df is not None:
                    df.to_csv(output.replace(".tsv", "") + "_enrichment_cluster%i.tsv" %cluster, sep="\t")
        if poscon:
            intersection_cluster_dict, common_genes = self.find_intersection_genes(poscon)
            for cluster in self.gene_dict:
                if len(intersection_cluster_dict[cluster]) >= 2:
                    print('intersection: cluster%i, len %i' %(cluster, len(intersection_cluster_dict[cluster])))
                    df = get_functional_enrichment(intersection_cluster_dict[cluster], cluster)
                    if df is not None:
                        df.to_csv(output.replace(".tsv", "") + "_enrichment_intersection_cluster%i.tsv" %cluster, sep="\t")

                if len(common_genes[cluster]) >= 2:
                    print('common: cluster%i, len %i' % (cluster, len(common_genes[cluster])))
                    df = get_functional_enrichment(common_genes[cluster], cluster)
                    if df is not None:
                        df.to_csv(output.replace(".tsv", "",) + "_enrichment_common_cluster%i.tsv" %cluster, sep='\t')




    def venn(self, posconsensus, output):
        poscon = Positional_Consensus(posconsensus)
        for cluster, mirna_list in self.seed_dict.items():
            cluster_genes_dict = {}
            if mirna_list:
                if len(mirna_list) >= 2 and len(mirna_list) <= 6:
                    for mirna in mirna_list:
                        genes = [x for x in poscon.mirna_dict.get(mirna, []) if x in self.gene_dict[cluster]]
                        cluster_genes_dict[mirna] = set(genes)

                    # Create a set of all genes
                    all_genes = set.union(*cluster_genes_dict.values())

                    # Create a binary matrix indicating gene presence for each miRNA
                    binary_matrix = pd.DataFrame(
                        {mirna: [1 if gene in cluster_genes_dict[mirna] else 0 for gene in all_genes] for mirna in
                         mirna_list}, index=list(all_genes))
                    venn_labels = {mirna: set(binary_matrix[binary_matrix[mirna] == 1].index) for mirna in mirna_list}

                    # Use venn function instead of venn3
                    venn.venn(venn_labels)

                    # if len(mirna_list) == 2:
                    #     # Create a Venn diagram
                    #     venn.venn2(venn_labels, names=mirna_list)
                    # elif len(mirna_list) == 3:
                    #     venn.venn3(venn_labels, names=mirna_list)
                    # elif len(mirna_list) == 4:
                    #     venn.venn4(venn_labels, names=mirna_list)
                    # elif len(mirna_list) == 5:
                    #     venn.venn5(venn_labels, names=mirna_list)
                    # elif len(mirna_list) == 6:
                    #     venn.venn6(venn_labels)
                    # else:
                    #     print("ERROR: mirna_list too large; cluster %i" %cluster)

                    # Set title
                    plt.title(f'Cluster {cluster} - Genes Targeted by miRNAs')

                    common_genes = set.intersection(*venn_labels.values())
                    contingency_table = [[len(venn_labels[mirna] - common_genes), len(common_genes)],
                                         [len(set(all_genes) - venn_labels[mirna]),
                                          len(all_genes) - len(venn_labels[mirna]) - len(common_genes)]]
                    _, p_value = fisher_exact(contingency_table)

                    # Display p-value on the plot
                    plt.text(0.5, -0.1, f'P-value: {p_value:.4f}', horizontalalignment='center',
                             verticalalignment='center', transform=plt.gca().transAxes)

                    # Save the plot to a file
                    plt.savefig(output.replace(".tsv", "") + "_vennplot_cluster%i.png" % cluster, dpi=300)

                    # Save the binary matrix to a CSV file
                    binary_matrix.to_csv(output.replace(".tsv", "") + "_vennplot_cluster%i.tsv" % cluster, sep='\t')

                    # Clear the current figure to avoid overlapping plots
                    plt.clf()
                else:
                    print("{} has {} mirnas".format(cluster, len(mirna_list)))



class Ematrix:
    def __init__(self,file):
        self.df = pd.read_csv(file, sep='\t')
        with open(file, "r") as fl:
            self.lines = fl.readlines()

    def get_median(self):
        numeric_columns = self.df.select_dtypes(include=['int', 'float']).apply(pd.to_numeric, errors='coerce')

        # Calculate the median along the columns
        median_values = numeric_columns.median(axis=1)

        # Add a new column to the DataFrame with the median values
        self.df['Median'] = median_values




class Positional_Consensus:
    def __init__(self, posconsensus):
        with open(posconsensus, 'r') as fl:
            self.lines = fl.readlines()
        self.mirna_dict = self.get_dict()

    def get_dict(self):
        mirna_dict = {}
        for line in self.lines:
            fields = line.strip('\n').split('\t')
            gene = fields[1]
            if "." in gene:
                gene = gene.split(".")[0]
            if "|" in gene:
                gene = gene.split("|")[0]
            mirna = fields[0]
            if mirna not in mirna_dict:
                mirna_dict[mirna] = [gene]
            else:
                mirna_dict[mirna].append(gene)
        return mirna_dict

class SeedMapping:
    def __init__(self, file):
        with open(file, 'r') as fl:
            self.lines = fl.readlines()
        self.seed_dict = self.get_dict()

    def get_dict(self):
        seed_dict = {}
        for line in self.lines:
            mirna = line.strip('\n').split('\t')[0]
            seed = line.strip('\n').split('\t')[1]
            if seed in seed_dict:
                seed_dict[seed].append(mirna)
            else:
                seed_dict[seed] = [mirna]

        return seed_dict

def find_common_genes(lists, min_occurrence=2):
    # Concatenate all lists into a single list
    all_genes = [gene for sublist in lists for gene in sublist]

    # Count occurrences of each gene
    gene_counts = Counter(all_genes)

    # Find genes that occur in at least min_occurrence lists
    common_genes = {gene for gene, count in gene_counts.items() if count >= min_occurrence}

    return common_genes

def seed_list_to_mirna(seeds, mirna_map_file):
    mapping = SeedMapping(mirna_map_file)
    mirna_list = []
    for seed in seeds:
        mirna_list = mirna_list + mapping.seed_dict[seed]

    return mirna_list


def get_functional_enrichment(gene_list, cluster=False):
    try:
        enrichment_df = stringdb.api.get_enrichment(gene_list, background_string_identifiers=None, species=9606,
                                               caller_identity='https://github.com/gpp-rnd/stringdb')
        return enrichment_df
    except:
        if cluster:
            print("cluster {} not retrieving enrichment".format(cluster))
            print(gene_list)
        pass




def network_import_directed(file, source, target):
    df = pd.read_csv(file, sep='\t')
    G = nx.from_pandas_edgelist(df, source=source, target=target, create_using=nx.DiGraph())
    return G

def network_stats(G, output):
    # Compute basic network statistics
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    density = nx.density(G)

    # Compute centrality measures
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)

    with open(output, 'wt') as out:
        out.write("NumNodes\tNumEdges\tDensity\tDegreeCentrality\tBetweennessCentrality\tClosenessCentrality\n")
        out.write(str(num_nodes) + "\t" + str(num_edges) + "\t" + str(density) + "\t" + str(degree_centrality) +
                  "\t" + str(betweenness_centrality) + "\t" + str(closeness_centrality) + "\n")


def network_community_analysis(G, output):
    # # Convert directed graph to undirected graph
    G_undirected = G.to_undirected()

    # Apply Louvain algorithm to the undirected graph
    partition = community.best_partition(G_undirected)

    print("Louvain Modularity:", community.modularity(partition, G_undirected))
    #
    # Visualize the Louvain communities
    colors = [partition[node] for node in G.nodes]
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
    plt.savefig(output + "_Louvain.png", dpi=300)

    # Perform Girvan-Newman community detection
    communities_generator = girvan_newman(G)
    top_level_communities = next(communities_generator)
    next_level_communities = next(communities_generator)

    print("Girvan-Newman Top Level Communities:", top_level_communities)
    print("Girvan-Newman Next Level Communities:", next_level_communities)

    # Visualize and save Girvan-Newman communities
    pos = nx.spring_layout(G)
    colors = [i for i, comm in enumerate(top_level_communities)]
    nx.draw(G, pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
    plt.savefig(output + "_Girvan-Newman_TL.png", dpi=300)
    plt.show()

    colors = [i for i, comm in enumerate(next_level_communities)]
    nx.draw(G, pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
    plt.savefig(output + "_Girvan-Newman_NL.png", dpi=300)
    plt.show()

    # Specify the value of k
    k_value = 10

    # Perform K-Clique Percolation community detection
    k_clique_communities = list(k_clique_communities(G, k_value))

    print(f"K-Clique Percolation Communities (k={k_value}):", k_clique_communities)

    # Visualize and save K-Clique Percolation communities
    pos = nx.spring_layout(G)
    for i, comm in enumerate(k_clique_communities):
        colors = [i] * G.subgraph(comm).number_of_nodes()
        nx.draw(G.subgraph(comm), pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
        plt.savefig(output + "_K-Clique.png", dpi=300)
        plt.show()

    # Perform Label Propagation community detection
    label_propagation_communities = list(label_propagation_communities(G))

    print("Label Propagation Communities:", label_propagation_communities)


    # Visualize and save Label Propagation communities
    pos = nx.spring_layout(G)
    for i, comm in enumerate(label_propagation_communities):
        colors = [i] * len(comm)
        nx.draw(G.subgraph(comm), pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
        plt.savefig(output + "_LabelPropagation.png", dpi=300)
        plt.show()

    # Perform Asynchronous Label Propagation community detection
    asyn_lpa_communities = list(asyn_lpa_communities(G))

    print("Asynchronous Label Propagation Communities:", asyn_lpa_communities)

    # Visualize and save Asynchronous Label Propagation communities
    pos = nx.spring_layout(G)
    for i, comm in enumerate(asyn_lpa_communities):
        colors = [i] * len(comm)
        nx.draw(G.subgraph(comm), pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
        plt.savefig(output + "_AsynLabelPropagation.png", dpi=300)
        plt.show()

    # Perform Greedy Modularity community detection
    greedy_modularity_communities = list(greedy_modularity_communities(G))

    print("Greedy Modularity Communities:", greedy_modularity_communities)

    # Visualize and save Greedy Modularity communities
    pos = nx.spring_layout(G)
    for i, comm in enumerate(greedy_modularity_communities):
        colors = [i] * len(comm)
        nx.draw(G.subgraph(comm), pos, node_color=colors, with_labels=True, cmap=plt.cm.Blues)
        plt.savefig(output + "_GreedyModularity.png", dpi=300)
        plt.show()


def network_kmeans(G, output, degree_centrality, num_clusters=3):
    degree_centrality_values = np.array(list(degree_centrality.values())).reshape(-1, 1)

    # Perform K-Means clustering
    kmeans = KMeans(n_clusters=num_clusters, random_state=42).fit(degree_centrality_values)
    cluster_labels = kmeans.labels_

    # Visualize the K-Means clusters
    pos = nx.spring_layout(G)
    nx.draw(G, pos, node_color=cluster_labels, with_labels=True, cmap=plt.cm.Blues)
    plt.savefig(output, dpi=300)

    # Example: Silhouette Score for K-Means clustering
    silhouette_avg = silhouette_score(degree_centrality_values, cluster_labels)
    print("Silhouette Score:", silhouette_avg)



