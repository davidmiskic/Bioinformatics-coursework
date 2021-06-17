import sys
import pickle
from os import path
import math

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
from scipy import spatial, cluster, stats
import matplotlib.pyplot as plt
from pandas import to_datetime, DataFrame
import seaborn

Entrez.email = "dm4018@fri.uni-lj.si"

"""
UTILITIES
"""
def jukes_cantor(p: float) -> float:
    """The Jukes-Cantor correction for estimating genetic distances.

    Parameters
    ----------
    p: float
        The proportional distance, i.e. the number of of mismatching symbols (Hamming
        distance) divided by the total sequence length.

    Returns
    -------
    float
        The corrected genetic distance.

    """
    return (-3 / 4) * math.log(1 - (4 * p / 3))

def load(organism_id: str) -> SeqRecord:
    """Load the NCBI record, use cached files if possible."""
    if not path.exists(path.join("data", f"{organism_id}.pkl.gz")):
        with Entrez.efetch(db="nucleotide", rettype="gb", id=organism_id) as handle:
            record = SeqIO.read(handle, "gb")
            with open(path.join("data", f"{organism_id}.pkl.gz"), "wb") as f:
                pickle.dump(record, f)
    else:
        with open(path.join("data", f"{organism_id}.pkl.gz"), "rb") as f:
            record = pickle.load(f)
    return record

def extractMetaData(identifier):
    record = load(identifier)
    label = []
    for feature in record.features:
        if feature.type == "source":
            date, host, country, organism = "/", ["Homo sapiens"], "/", "/"
            if "collection_date" in feature.qualifiers.keys(): date = feature.qualifiers["collection_date"]
            if "host" in feature.qualifiers.keys(): host = feature.qualifiers["host"]
            if "country" in feature.qualifiers.keys(): country = feature.qualifiers["country"]
            if "organism" in feature.qualifiers.keys(): organism = feature.qualifiers["organism"]
            label = [organism[0], host[0], date[0], country[0]]

    return label

def findMin(distances):
    rows, cols = distances.shape
    minel = 10e6
    mini, minj = 0,0
    for i in range(0, rows):
        for j in range(0, cols):
            if i != j and i > j and distances[i][j] < minel:
                mini, minj = i, j
                minel = distances[i][j]
    return (minel, mini, minj)

def UPGMA(distances):
    """Unweighted pair group method with arithmetic mean (UPGMA) agglomerative clustering.
    
    Parameters
    ----------
    distances: np.ndarray
        A two dimensional, square, symmetric matrix containing distances between data
        points. The diagonal is zeros.
        
    Returns
    -------
    np.ndarray
        The linkage matrix, as specified in scipy. Briefly, this should be a 2d matrix
        each row containing 4 elements. The first and second element should denote the
        cluster IDs being merged, the third element should be the distance, and the
        fourth element should be the number of elements within this new cluster. Any
        new cluster should be assigned an incrementing ID, e.g. after the first step
        where the first two points are merged into a cluster, the new cluster ID should
        be N, then N+1, N+2, ... in subsequent steps.
    
    """
    rows, cols = distances.shape
    # last used cluster ID, new ones should be N+1...
    lastCluster = rows

    # key is merged cluster id's, value is number of elements
    clusterIDs = {}
    composition = {}
    for i in range(0, rows):
        clusterIDs[str(i)] = 1
        composition[i] = [i]

    linkage = []
    distancesOriginal = np.copy(distances)
    maxel = distancesOriginal.max()
    labels = [i for i in range(rows)]
    while rows > 1:
        # print(labels)
        # print(distances)

        # mini je VRSTICA, minj je STOLPEC
        minel, mini, minj = findMin(distances)
        # print(mini, minj)
        # axis 0=row, 1=col
        # new distances are average of joined
        newdists = [0 for i in range(0, rows)]

        for j in range(0, rows): newdists[j] = (distances[mini][j] * len(composition[labels[mini]]))
        for i in range(0, cols): newdists[i] = (newdists[i]+distances[minj][i]*len(composition[labels[minj]])) / (len(composition[labels[mini]]) + len(composition[labels[minj]]))
        # correct diagonal to 0
        newdists[mini] = 0
        newdists[minj] = 0

        newcls = [min(labels[mini], labels[minj]) , max(labels[mini], labels[minj]),
                  minel, clusterIDs[str(labels[mini])] + clusterIDs[str(labels[minj])]]
        clusterIDs[str(lastCluster)] = newcls[3]
        linkage.append(newcls)
        # delete ith row and column and change jth col values
        for j in range(0, len(newdists)):
            distances[j][minj] = newdists[j]
            distances[minj][j] = newdists[j]
            # distances[j][mini] = maxel + 1
            # distances[mini][j] = maxel + 1
        #distances = np.append(distances, newdists, 1)

        distances = np.delete(np.delete(distances, mini, 0), mini, 1)
        composition[lastCluster] = composition[labels[mini]] + composition[labels[minj]]
        labels.pop(mini)
        labels[minj] = lastCluster
        lastCluster += 1
        rows, cols = distances.shape

        # print(linkage)
        # print("\n")
    # print(clusterIDs)
    # print(composition)
    return np.array(linkage)

def dist(s1, s2):
    return jukes_cantor(sum([1 for i in range(0, len(s1)) if s1[i] != s2[i]])/(float(len(s1))))

"""
PROCESSING
"""

def parseDataToDistMatrix(path="data/sars-cov-2.fa"):
    # sequences already aligned - all of the same length
    records = list(SeqIO.parse(path, "fasta"))
    distances = np.zeros((len(records), len(records)), np.float)
    labels = []

    for i in range(0, len(records)):
        temp = extractMetaData(records[i].id)
        labels.append(str(records[i].id) + ", " + temp[-1] + ", " + temp[2])
    print(labels)
    # print(records[107].id)
    # extractMetaData(records[107].id)
    for i in range(0, len(records)):
        print(i)
        for j in range(i, len(records)):
            if i != j: distances[i][j] = dist(records[i].seq, records[j].seq)
    return (distances, labels)


def dendrogram(distances, lab):
    linkage = UPGMA(distances + distances.transpose())
    plt.figure(figsize=(200, 400))
    fig, axes = plt.subplots()
    cluster.hierarchy.dendrogram(linkage, orientation="left",  labels=lab
                                 , color_threshold=550, leaf_font_size=2.5)
    plt.savefig("problem2.png", bbox_inches = 'tight', dpi=500)
    print("dendrogram complete")


def mutationRate(path="data/sars-cov-2.fa", name="problem3a.png", zoo=False):
    records = list(SeqIO.parse(path, "fasta"))
    collectionDates = []
    distances = []
    dateDistances = []
    for i in range(0, len(records)):
        temp = extractMetaData(records[i].id)
        collectionDates.append(to_datetime(temp[2]))
    earliest = collectionDates.index(min(collectionDates))
    for i in range(0, len(records)):
        distances.append(dist(records[earliest].seq, records[i].seq))
        dateDistances.append((collectionDates[i] - collectionDates[earliest]).days)

    xs = [] # days difference
    ys = [] # distance
    collectionDatesForZoo = collectionDates.copy()
    collectionDates.pop(earliest)
    dateDistances.pop(earliest)
    distances.pop(earliest)
    xs.append(0)
    ys.append(0)
    while len(collectionDates) != 0:
        next = collectionDates.index(min(collectionDates))
        xs.append(dateDistances[next])
        ys.append(distances[next])
        collectionDates.pop(next)
        dateDistances.pop(next)
        distances.pop(next)

    if name == "problem3bHIV.png":
        xs = xs[1:]
        ys = ys[1:]

    model = stats.linregress(xs, ys)
    print(model)
    modelYs = [model.slope*xs[i] + model.intercept for i in range(0, len(xs))]
    plt.figure(figsize=(20, 20))
    plt.plot(xs, ys, ".")
    plt.plot(xs, modelYs, "-")
    plt.xlabel("distance = " + str(model.slope) + "*dayDifference " + "+ " + str(model.intercept))
    plt.ylabel("JC-Hamming distance")
    plt.title("Linear regression for days X distance")

    if zoo:
        recordsZoo = list(SeqIO.parse("data/sars-cov-2-animals.fa", "fasta"))
        distZoo = [dist(records[earliest].seq, recordsZoo[i].seq) for i in range(0, len(recordsZoo))]
        dateDistZoo = []
        labelsZoo = []
        for i in range(0, len(recordsZoo)):
            temp = extractMetaData(recordsZoo[i].id)
            labelsZoo.append(temp[1])
            dateDistZoo.append((to_datetime(temp[2]) - (collectionDatesForZoo[earliest])).days)
        plt.scatter(dateDistZoo, distZoo, c="r")
        for x,y in zip(dateDistZoo, distZoo): plt.annotate(labelsZoo.pop(0), (x, y), textcoords="offset points", xytext=(0, 10))

    plt.savefig(name)



"""
ORF FINDER
"""
def isNotNested(start, end, others):
    for orf in others:
        if orf[0] <= start and orf[1] >= end: return False
    return True


def translate_to_protein(seq):
    if len(seq) % 3 != 0: seq = seq[0:len(seq)-(len(seq) % 3)]
    standardTable = { "TTC":"F", "TTT":"F",
                      "CTC": "L", "CTG": "L", "CTT": "L", "CTA":"L", "TTA":"L", "TTG":"L",
                      "ATA":"I", "ATC":"I", "ATT":"I",
                      "ATG":"M",
                      "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
                      "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "AGC":"S", "AGT":"S",
                      "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                      "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
                      "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                      "TAC": "Y", "TAT": "Y",
                      "CAC": "H", "CAT": "H",
                      "CAA": "Q", "CAG": "Q",
                      "AAC":"N", "AAT":"N",
                      "AAA":"K", "AAG":"K",
                      "GAC": "D", "GAT": "D",
                      "GAA": "E", "GAG": "E",
                      "TGC": "C", "TGT": "C",
                      "TGG": "W",
                      "AGA":"R", "AGG":"R",
                      "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
                      "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G"
                      #"TAA":"*", "TAG":"*", "TGA":"*"
    }
    translation = ""
    for i in range(0, len(seq) - 2, 3):
        try: translation += standardTable[seq[i:i+3]]
        except KeyError: pass

    return translation


def find_orfs(sequence, start_codons = ["ATG"], stop_codons = ["TGA"]):
    """Find possible ORF candidates in a single reading frame.
    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)
    """
    found = []
    for i in range(0, len(sequence) - 2, 3):
        if sequence[i:i + 3] in start_codons:
            for j in range(i, len(sequence) - 2, 3):
                if sequence[j:j + 3] in stop_codons:
                    # no nested orfs
                    if isNotNested(i, j+3, found): found.append((i,j+3))
                    break
    return found

def find_all_orfs(sequence, start_codons = ["ATG"], stop_codons = ["TGA"]):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.
    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.
    """
    # TODO
    found = []
    revCom = sequence.reverse_complement()
    for i in range(0, 3, 1):
        orfs = find_orfs(sequence[i:], start_codons, stop_codons)
        for orf in orfs: found.append((1, orf[0]+i, orf[1]+i))
    for i in range(0, 3, 1):
        orfs = find_orfs(revCom[i:], start_codons, stop_codons)
        for orf in orfs: found.append((-1, len(sequence) - (orf[1] + i), len(sequence) - (orf[0] + i)))
    return found


def snps(path="data/sars-cov-2.fa"):
    records = list(SeqIO.parse(path, "fasta"))
    # distribution of nucleotides across all sequences and their entropy

    dists = []
    l = len(records[0].seq)
    ls = len(records)
    for i in range(0, l):
        aa = sum([records[a].seq[i] == "A" for a in range(0, len(records))]) / float(ls)
        cs = sum([records[a].seq[i] == "C" for a in range(0, len(records))]) / float(ls)
        ts = sum([records[a].seq[i] == "T" for a in range(0, len(records))]) / float(ls)
        gs = sum([records[a].seq[i] == "G" for a in range(0, len(records))]) / float(ls)
        dists.append([aa, cs, ts, gs])
    ents = [stats.entropy(dist, base=2) for dist in dists]

    # f = open("entropies.csv", "w")
    # cnt = 1
    # for ent in ents:
    #     f.write(str(cnt) + "," + str(ent) + "\n")
    #     cnt += 1
    # f.close
    # print("file written")
    xs = []
    ys = []
    for i in range(0, len(ents)):
        if ents[i] > 0:
            ys.append(ents[i])
            xs.append(i+1)

    # plot entropies
    #ax.bar(x=xs, height=ys)
    # for i in range(0, len(xs)):
    #     plt.text(xs[i], ys[i], str(round(ys[i], 2)), ha="center")
    plt.figure(figsize=(150, 10))
    #ax = seaborn.barplot(y="Location" ,x="Entropy", data=DataFrame({"Entropy":ys, "Location":xs}), orient="h")
    ax = seaborn.barplot(x="Location", y="Entropy", data=DataFrame({"Entropy":ys, "Location":xs}))

    for p in ax.patches:
        ax.annotate(round(p.get_height(), 2),
           (p.get_x() + p.get_width() / 2., p.get_height()),
           ha="center", va="center", xytext=(0, 9), textcoords="offset points")

    plt.xticks(rotation=90)
    plt.savefig("mpl.png", bbox_inches = 'tight')
    #plt.show()



def snpsBetweenSequences():
    records = list(SeqIO.parse(path, "fasta"))
    # compare longest ORF on each sequence to reference
    f = open("sequence.txt", "r")
    f.readline()
    firstorf1ab1 = f.readline()

    longestORF = (1, 172, 13963)  # they all are the same in this case
    ORFtranslation = []
    for rec in records:
        print(rec.id)
        # found = find_all_orfs(rec.seq)
        # maxf = found[0]
        # for f in found:
        #    if f[2] - f[1] > maxf[2] - maxf[1]: maxf = f
        # ORFtranslation.append(translate_to_protein( rec.seq[longestORF[1]:longestORF[2]]) )
        ORFtranslation.append(translate_to_protein(rec.seq[longestORF[1]:longestORF[2]]))
        # longestORF = maxf
    referencetranslation = firstorf1ab1[0:len(ORFtranslation[0])]

    diffmatrix = []
    for prot in ORFtranslation:
        protdiff = ""
        for i in range(0, len(prot)):
            if prot[i] == referencetranslation[i]:
                protdiff += "-"
            else:
                protdiff += prot[i]
        diffmatrix.append(protdiff)

    # f = open("ab1comparison.txt", "w")
    # f.write(referencetranslation + "\n")
    # for i in range(0, len(records)):
    #     f.write(">" + records[i].id + "\n")
    #     f.write(diffmatrix[i] + "\n")
    # f.close()
    return diffmatrix

if __name__ == '__main__':
    phylo = False
    mutation = False
    mutationOthers = False
    mutationZoo = False
    observeSNP = False


    if phylo:
        dis, lab = parseDataToDistMatrix()
        dendrogram(dis, lab)
    if mutation:
        mutationRate()
    if mutationOthers:
        mutationRate(path="data/ebola.fa", name="problem3bEbola.png")
        mutationRate(path="data/hiv.fa", name="problem3bHIV.png")
    if mutationZoo:
        mutationRate(name="problem4.png", zoo=True)
    if observeSNP:
        snps()
        snpsBetweenSequences()

"""
TESTIRANJE

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage
print(linkage(squareform(sample.transpose() + sample), "average"))

sample = np.array([
    [0, 19, 27, 8, 33, 18, 13],
    [0, 0, 31, 18, 36, 1, 13],
    [0, 0, 0, 26, 41, 32, 29],
    [0, 0, 0, 0, 31, 17, 14],
    [0, 0, 0, 0, 0, 35, 28],
    [0, 0, 0, 0, 0, 0, 12],
    [0, 0, 0, 0, 0, 0, 0]
], np.float)
a = (UPGMA(sample.transpose() + sample))
plt.figure(figsize=(20,10))
fig, axes = plt.subplots()
cluster.hierarchy.dendrogram(a, orientation="left"
                                          ,color_threshold=550, leaf_font_size=8)
plt.show()

"""


