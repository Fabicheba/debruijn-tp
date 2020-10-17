""" Script qui permet d'assembler le génome de de l'entérovirus,
génome très court (7408 nucléotides), linéaire et non segmenté..
Nous impléménterons donc un graphe de de Bruijn
----------------------------------------------------------------
Etudiante: BARRY Fatoumata Binta M2 BI Université Paris Diderot
UE: Genomique / Assemblage
Année : 2020-2021
----------------------------------------------------------------

"""

# Modules utiles
import argparse
import networkx as net
import matplotlib.pyplot as plt
# Fonctions

def argument():
    """ fonctions qui permet de d'ajouter des arguments
    pour lancer le script debruijn.
    ----------
    Returns
    -------
    tuple
        fastq: fichier d'entree
        k: taille des kmer (entier positif)
        contig: fichier output avec les contigs
    """

    #Declaration des arguments
    arg = argparse.ArgumentParser()
    arg.add_argument("-i", "--fastq",  help="fichier fastq single end.", type=str, required=True)
    arg.add_argument("-k", "--kmer", help="taille des kmers", type=int, required=True, default=21)
    arg.add_argument("-o", "--contig", help="fichier output avec les contigs, nom à choisir",
                        type=str, required=True)

    args = arg.parse_args()
    fastq_file = args.fastq
    kmer_num = args.kmer
    contig_file = args.contig

    return fastq_file, kmer_num, contig_file


# Création du graphe de de Bruijn
# -----------------------------------
# a. Identification des k-mer unique

def read_fastq(fastq_file):
    """ fonctions qui prend en argument un fichier fastq
    et return un generateur de sequence

    Parametre
    ----------
    fastq_file: fichier fastq donné en argument au script
    Returns
    -------
    generateur de sequences: str

    """
    with open(fastq_file, "r") as fast:
        for ligne in fast:
            if not ligne.startswith("@") and (not ligne.startswith("+")) and (not ligne.startswith("J")):
                yield ligne[:-1]


def cut_kmer(seq, k):
    """" fonctions qui prend en argument une sequence, une taille
            de kmer et qui returne un generateur de kmer

    Parametres
    ----------
    seq: sequence issue du fichier fastq
    k: int, taille de kmer
    Returns
    -------
    generateur de kmer: str

    """
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]


def build_kmer_dict(fastq_file, k):
    """" fonction qui prend en argument un fichier fastq, une taille
        de kmer et et retourne un dictionnaire ayant pour clé le k-mer
        et pour valeur le nombre d’occurrence de ce k-mer

    Parametres
    ----------
    fasta_file: fichier fastq
    k: int, taille de kmer
    Returns
    -------
    dict_kmer: dict

    """
    dict_kmer = {}

    for seq in read_fastq(fastq_file): # appel du generateur de sequence
        for kmer in cut_kmer (seq, k): # appel du generateur de kmer
            if kmer not in dict_kmer:
                dict_kmer[kmer] =1
            else:
                dict_kmer[kmer] += 1

    return dict_kmer


# b. Construction de l’arbre de de Bruijn

def build_graph(dict_kmer):
    """" fonction qui prend en entrée un dictionnaire de k-mer
         et crée l’arbre de k-mers préfixes et suffixes décrit
         dans la fonction build_kmer_dict.

    Parametres
    ----------
    dict_kmer: dictionnaire d'occurence des kmers
    -------
    graphe: dico avec comme clé un tuple (prefixe, suffixe)
                et comme valeur le poids(occurence du kmer)

    """
    graphe=net.DiGraph()
    for i in dict_kmer:
        prefixe = i[:-1]
        suffixe = i[1:]
        poids = dict_kmer[i]
        graphe.add_edge(prefixe, suffixe, weight= poids)
    net.draw(graphe)
    plt.show()
    return graphe


# 2 Parcours du graphe de de Bruijn

"""def get_starting_nodes(graphe):
    pref_suff=list(graphes.keys())
    entree=[]
    for i in len(range(pref_suf)-1):
        for j in len(range(1,pref_suff)):
            if pref_suff[i][0] != pref_suff[j][0]"""

# Programme principal
if __name__ == "__main__":
    fastq, k, contig = argument()
    kmer = build_kmer_dict(fastq, k)
    print(kmer)
    print(build_graph(kmer))
