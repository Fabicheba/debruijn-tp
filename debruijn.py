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
import os
from statistics import stdev

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
    graphe: graphe avec le module Networkx(prefixe, suffixe)
                et comme poids l'occurence du kmer

    """
    graphe=net.DiGraph()
    for i in dict_kmer:
        prefixe = i[:-1]
        suffixe = i[1:]
        poids = dict_kmer[i]
        graphe.add_edge(prefixe, suffixe, weight= poids)
    net.draw(graphe)
    #plt.show()
    return graphe


# 2 Parcours du graphe de de Bruijn

def get_starting_nodes(graphe):
    """" fonction qui prend en entrée un graphe et 
        retourne une liste de noeuds d'entrée
    Parametres
    ----------
    graphe: graph de Networkx
    -------
    entree: liste des noeuds d'entree

    """
    noeuds_entree = []
    for noeud in graphe.nodes:
        if len(list(graphe.predecessors(noeud))) == 0:
            noeuds_entree.append(noeud)
    return noeuds_entree



def get_sink_nodes(graphe):
    """" fonction qui prend en entrée un graphe et 
        retourne une liste de noeuds de sortie
    Parametres
    ----------
    graphe: graph de Networkx
    -------
    sortie: liste des noeuds de sortie

    """
    noeuds_sortie = []
    for noeud in graphe.nodes:
        if len(list(graphe.successors(noeud))) == 0:
            noeuds_sortie.append(noeud)
    return noeuds_sortie


def get_contigs(graphe, noeuds_entree, noeuds_sortie):
    """" fonction qui prend un graphe, une liste de noeuds
        d’entrée et une liste de sortie et retourne une 
        liste de tuple(contig, taille du contig)
    Parametres
    ----------
    graphe: graph de Networkx
    noeuds_entree: liste des noeuds d'entree
    noeuds_sortie: liste des noeuds de sortie
    Return
    -------
    contigs: liste de tuple (contig, taille du contig)

    """
    liste_contigs = []
    for entree in noeuds_entree:
        for sortie in noeuds_sortie:
            for chemin in net.all_simple_paths(graphe, source=entree, target=sortie):
                contig = chemin[0]
                for i in range(1, len(chemin)):
                    contig += chemin[i][-1]
                liste_contigs.append((contig, len(contig)))
    return liste_contigs
        

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(liste_contigs, contigs_file):
    with open(contigs_file, "w") as f_contig:
        num_contig = 0
        for contig in liste_contigs:
            num_contig += 1
            f_contig.write{">contig_{}  len = {}\n".format(num_contig, contig[1])}
            f_contig.write(fill(contig[0]) + "\n")


# 3. Simplification du graphe de de Bruijn
# a. Résolution des bulles

def std(val_liste):
    """" fonction qui prend une liste de valeurs et
         qui retourne l'ecart type
    Parametres
    ----------
    val_list: liste de valeurs
    Return
    -------
    ecart_type: ecart type de la valeur de la liste en entree
    
    """
    return stdev(val_liste)


def path_average_weigth(graphe, chemin):
    """" fonction qui prend un graphe et un chemin et
        retourne un poids moyen
    Parametres
    ----------
    graphe: Graphe Networkx, issu de la fonction build_graphe
    chemin: chemin dans le graphe
    Return
    -------
    poids_moyen: poids moyen de tous les chemins

    """
    poids = 0
    for i in range(len(chemin)-1):
        poids += graphe.edges[chemin[i], chemin[i+1]]["weigth"]
    return poids/(len(path)-1) # poids moyen


def remove_paths(graphe, liste_chemin, delete_entry_node,
                    delete_sink_node):
    """" fonction qui permet de nettoyer le graphe,
        en supprimant les chemins indésirables
    Parametres
    ----------
    graphe: Graphe Networkx, issu de la fonction build_graphe
    liste_chemin: chemin dans le graphe
    delete_entry_node: booleen
        pour indiquer si les noeuds d’entrée seront supprimés
    delete_sink_node: booleen
        pour indiquer si les noeuds
         de sortie seront supprimés
    Return
    -------
    graphe: nettoyé des chemins indésirables

    """
    entry = sink = 1
    if delete_entry_node:
        entry = 0
    if delete_sink_node:
        sink = 0
    for chemin in liste_chemin:
        for i in range(entry, len(chemin) - sink):
            graphe.remove_node(chemin[i])
    return graphe


# Programme principal
if __name__ == "__main__":
    fastq, k, contig = argument()
    kmer = build_kmer_dict(fastq, k)
    #print(kmer)
    graphe= build_graph(kmer)
    entree=get_starting_nodes(graphe)
    sortie= get_sink_nodes(graphe)
    contig = get_contigs(graphe, entree, sortie)
    print(entree, sortie)
    print(contig)