import pickle
from os import path
from typing import Tuple, Generator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial, cluster

import multiprocessing as mp

"""
GENERAL
"""
Entrez.email = "dm4018@fri.uni-lj.si"

def hamming(s1, s2):
    return sum([1 for i in range(0, len(s1)) if s1[i] != s2[i]])

def blastScore(s1, s2):
    return sum([1 for i in range(0, len(s1)) if s1[i] == s2[i]])

def scoring(x, y):
    if [x, y] == ["*", "*"]: return 1
    elif "*" in [x, y]: return -4
    else:
        try:
            return blosum62[(x, y)]
        except KeyError:
            return blosum62[(y, x)]

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

def computeAlignments_writeOut(extraction, scoringFunction):
    scorematrix = np.zeros((21, 21), np.int)
    allalignments = {}
    f = open("alignments.txt", "w")
    for i in range(len(extraction.keys())):
        allalignments[i] = {}
        for j in range(0, len(extraction.keys())):
            print(i, j)
            alignmentSeq1, alignmentSeq2, score = global_alignment(str(extraction[list(extraction.keys())[i]]),
                                                                   str(extraction[list(extraction.keys())[j]]), scoringFunction)
            allalignments[i][j] = [alignmentSeq1, alignmentSeq2, score]
            scorematrix[i][j] = score
            f.write(">" + str(list(extraction.keys())[i]) + "<->" + str(list(extraction.keys())[j]) + "\n")
            f.write(alignmentSeq1 + "\n" + alignmentSeq2 + "\n" + str(score) + "\n")
            f.flush()
    f.close()
    np.savetxt("scoreMatrix.txt", scorematrix, fmt="%.0f", delimiter=",", header=str(list(extraction.keys())))
    return scorematrix

def loadAlignments(file):
    allalignments = {}
    f = open(file, "r")
    l = f.readline()
    while len(l) > 0:
        if l[0] == ">":
            i, j = l[1:].split("<->")
            j = j.replace("\n", "")
            if i not in allalignments.keys():  allalignments[i] = {}
            allalignments[i][j] = [f.readline().replace("\n", ""), f.readline().replace("\n", ""), f.readline().replace("\n", "")]
        l = f.readline()
    f.close()
    return allalignments

def blastQuality(s1, s2):
    # coverage of second with first
    return sum(
        [1 for i in range(0, len(s1)) if s1[i] == s2[i] and [s1] != "-"]
    ) /len(s2.replace("-", ""))

"""
ALIGNMENT
"""
def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    penalty = -1
    scorematrixH = np.zeros((len(seq1)+1, len(seq2)+1))
    highest = [0,0]
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            align = scorematrixH[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            gapfirst = scorematrixH[i][j-1] + penalty
            gapsecond = scorematrixH[i-1][j] + penalty
            scorematrixH[i][j] = max(0, align, gapfirst, gapsecond)
            if scorematrixH[i][j] > scorematrixH[highest[0]][highest[1]]:
                highest = [i, j]
    alignment1 = ""
    alignment2 = ""
    i, j = highest
    while scorematrixH[i][j] != 0:
        # diagonal: match, up: gap in first seq, left: gap in second seq
        diag = scorematrixH[i-1][j-1]
        up = scorematrixH[i-1][j]
        left = scorematrixH[i][j-1]
        if diag == max(diag, up, left):
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i = i - 1
            j = j - 1
        elif left == max(diag, up, left):
            alignment1 = "*" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j = j - 1
        else:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "*" + alignment2
            i = i - 1
    # print(seq1)
    # print(alignment1)
    # print(alignment2)
    # print(seq2)
    score = sum([scoring_function(alignment1[i], alignment2[i]) for i in range(0, len(alignment1))])
    # print(score)
    return alignment1.replace("*", "-"), alignment2.replace("*", "-"), score


def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    penalty = -1
    alignmentMatrix = np.zeros((len(seq1)+1, len(seq2)+1))
    for i in range(0, len(seq1)+1): alignmentMatrix[i][0] = penalty*i
    for j in range(0, len(seq2)+1): alignmentMatrix[0][j] = penalty*j

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            match = alignmentMatrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            delete = alignmentMatrix[i-1][j] + penalty
            insert = alignmentMatrix[i][j-1] + penalty
            alignmentMatrix[i][j] = max(match, delete, insert)

    alignmentSeq1 = ""
    alignmentSeq2 = ""
    i = len(seq1)
    j = len(seq2)
    while i > 0 and j > 0:
        if i > 0 and j > 0 and alignmentMatrix[i,j] == alignmentMatrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]):
            alignmentSeq1 = seq1[i-1] + alignmentSeq1
            alignmentSeq2 = seq2[j-1] + alignmentSeq2
            i = i - 1
            j = j - 1
        elif i > 0 and alignmentMatrix[i][j] == alignmentMatrix[i-1][j] + penalty:
            alignmentSeq1 = seq1[i-1] + alignmentSeq1
            alignmentSeq2 = "*" + alignmentSeq2
            i = i - 1
        elif j > 0 and alignmentMatrix[i][j] == alignmentMatrix[i][j-1] + penalty:
            alignmentSeq1 = "*" + alignmentSeq1
            alignmentSeq2 = seq2[j-1] + alignmentSeq2
            j = j - 1
    # print(i, j)
    # print(seq1)
    # print(alignmentSeq1)
    # print(alignmentSeq2)
    # print(seq2)
    score = sum([scoring_function(alignmentSeq1[i], alignmentSeq2[i]) for i in range(0, len(alignmentSeq1))])
    #print(score)
    return alignmentSeq1.replace("*", "-"), alignmentSeq2.replace("*", "-"), score

"""
PROCESSING
"""
def extractSpikeProtein(identifier):
    record = load(identifier)
    relevant = []
    foundLocations = set()
    for feature in record.features:
        if feature.type == "CDS":
            if ("gene" in feature.qualifiers.keys() and "S" in feature.qualifiers["gene"]) or "spike" in feature.qualifiers["product"][0]:
                if "translation" not in feature.qualifiers.keys(): print("warning! No translation for CDS.")
                else:
                    relevant.append(feature)
                    foundLocations.add((feature.location.start, feature.location.end))
    for feature in record.features:
        if feature.type == "gene":
            if "gene" in feature.qualifiers.keys() and "S" in feature.qualifiers["gene"]:
                if (feature.location.start, feature.location.end) not in foundLocations:
                    if feature.location.strand == 1: geneSequence = record.seq[feature.location.start:feature.location.end]
                    else: geneSequence = record.seq[feature.location.start:feature.location.end].reverseComplement()
                    feature.mytranslation = geneSequence.translate()
                    relevant.append(feature)
    result = []
    for feature in relevant:
        if feature.type == "CDS": result.append(feature.qualifiers["translation"][0])
        else: result.append(feature.mytranslation)
    return result


def extractCDS(identifier, organism):
    record = load(identifier)
    relevant = []
    for feature in record.features:
        if feature.type == "CDS":
            seq = ""
            name = organism
            if "translation" not in feature.qualifiers.keys(): print("warning! No translation for CDS.")
            else: seq = feature.qualifiers["translation"][0]
            if "gene" in feature.qualifiers.keys(): name = name + " " + feature.qualifiers["gene"][0]
            if "product" in feature.qualifiers.keys(): name = name + " " + feature.qualifiers["product"][0]
            relevant.append([name, seq])
    return relevant


def dendrogram(alignments, file, distance = "H"):
    distanceMatrix = np.zeros((len(alignments),len(alignments)))
    for i in range(len(alignments.keys())):
        for j in range(i+1, len(alignments.keys())):
            namei = list(alignments.keys())[i]
            namej = list(alignments.keys())[j]
            if distance == "H": distanceMatrix[i][j] = hamming(alignments[namei][namej][0], alignments[namei][namej][1])
            else: distanceMatrix[i][j] = len(alignments[namei][namej][0]) - blastScore(alignments[namei][namej][0], alignments[namei][namej][1])
            distanceMatrix[j][i] = distanceMatrix[i][j]
    # for i in distanceMatrix: print(list(i))
    # print(np.allclose(distanceMatrix, distanceMatrix.T))
    np.savetxt("denddistancematrix.txt", distanceMatrix, fmt="%.0f", delimiter=",",
               header=str(list(alignments.keys())))
    distanceMatrix = spatial.distance.squareform(distanceMatrix)
    clusters = cluster.hierarchy.linkage(distanceMatrix, method="average")

    plt.figure(figsize=(20,10))
    fig, axes = plt.subplots()
    dendro = cluster.hierarchy.dendrogram(clusters, labels=list(alignments.keys()), orientation="left"
                                          ,color_threshold=550, leaf_font_size=8)
    #plt.show()
    plt.savefig(file, bbox_inches = 'tight', dpi=500)


def parallelMatrixGlobal(i, j, extraction):
    print(i, j)
    alignmentSeq1, alignmentSeq2, score = global_alignment(str(extraction[list(extraction.keys())[i]]),
                                                           str(extraction[list(extraction.keys())[j]]), scoring)
    return [(i, j, score), list(extraction.keys())[i], list(extraction.keys())[j], alignmentSeq1, alignmentSeq2]


def parallelMatrixLocal(i, j, extraction):
    print(i, j)
    alignmentSeq1, alignmentSeq2, score = local_alignment(str(extraction[list(extraction.keys())[i]]),
                                                           str(extraction[list(extraction.keys())[j]]), scoring)
    return [(i, j, score), list(extraction.keys())[i], list(extraction.keys())[j], alignmentSeq1, alignmentSeq2]

# from Bio.SubsMat import MatrixInfo as matlist
# matrix = matlist.blosum62
# print(matrix.keys())
#
# from Bio.Align import substitution_matrices
# matrix2 = substitution_matrices.load("blosum62")
# print(matrix2)

def computeLocal(extraction):
    f = open("alignmentsLocal.txt", "w")
    dim = len(extraction)
    pool = mp.Pool(mp.cpu_count())
    parallelresult = pool.starmap(parallelMatrixLocal, [(i, j, extraction) for i in range(0, dim) for j in range(i, dim)])
    pool.close()  # no more tasks
    pool.join()  # wait for the remaining tasks to complete

    scorematrix = np.zeros((dim, dim), np.int)
    for res in parallelresult:
        i, j, score = res[0]
        scorematrix[i][j] = score
        f.write(">" + res[1] + "<->" + res[2] + "\n")
        f.write(res[3] + "\n" + res[4] + "\n" + str(score) + "\n")
        f.flush()
    f.close()
    np.savetxt("scoreMatrixLocal.txt", scorematrix, fmt="%.0f", delimiter=",", header=str(list(extraction.keys())))

def computeGlobal(extraction):
    f = open("alignments.txt", "w")
    pool = mp.Pool(mp.cpu_count())
    dim = len(extraction)
    parallelresult = pool.starmap(parallelMatrixGlobal,
                                  [(i, j, extraction) for i in range(0, dim) for j in range(i, dim)])

    # [(i, j, score), list(extraction.keys())[i], list(extraction.keys())[j], alignmentSeq1, alignmentSeq2]

    pool.close()  # no more tasks
    pool.join()  # wait for the remaining tasks to complete
    scorematrix = np.zeros((dim, dim), np.int)
    for res in parallelresult:
        i, j, score = res[0]
        scorematrix[i][j] = score
        f.write(">" + res[1] + "<->" + res[2] + "\n")
        f.write(res[3] + "\n" + res[4] + "\n" + str(score) + "\n")
        f.flush()
    f.close()
    np.savetxt("scoreMatrix.txt", scorematrix, fmt="%.0f", delimiter=",", header=str(list(extraction.keys())))

def miniBLAST(reference_genomes, query, orf_candidates, file=""):
    extraction2 = {}
    for name in reference_genomes:
        temp = extractCDS(accession_codes[name], name)
        for item in temp:
            #
            extraction2[item[0]] = item[1]

    queryGenome = load(accession_codes[query]).seq
    queryORFS = {}
    for orf in orf_candidates:
        s, e = orf_candidates[orf][1], orf_candidates[orf][2]
        queryORFS[orf] = queryGenome[s:e]
        if orf_candidates[orf][0] == -1: queryORFS[orf].reverse_complement()
        queryORFS[orf] = str(queryORFS[orf].translate())
        extraction2[orf] = queryORFS[orf]

    pool = mp.Pool(mp.cpu_count())
    parallelresult = pool.starmap(parallelMatrixLocal, [(i, j, extraction2) for i in range(0, len(extraction2)) for j in
                                                        range(i, len(extraction2))])
    pool.close()  # no more tasks
    pool.join()  # wait for the remaining tasks to complete

    scorematrix = np.zeros((len(extraction2), len(extraction2)), np.int)
    f = open("miniBlast" + f + ".txt", "w")
    for res in parallelresult:
        i, j, score = res[0]
        scorematrix[i][j] = blastScore(res[3], res[4])
        f.write(">" + res[1] + "<->" + res[2] + "\n")
        f.write(res[3] + "\n" + res[4] + "\n" + str(scorematrix[i][j]) + "\n")
        f.flush()
    f.close()
    np.savetxt("scoreMatrixMiniBlast" + file + ".txt", scorematrix, fmt="%.0f", delimiter=",", header=str(list(extraction2.keys())))

def processBLASTResults(allalignments, orf_candidates):
    scoresblast = {}
    for k1 in allalignments.keys():
        if k1 in orf_candidates: continue
        for k2 in allalignments[k1].keys():
            if k2 in orf_candidates:
                if len(allalignments[k1][k2][0]) != len(allalignments[k1][k2][1]):
                    raise Warning
                else:
                    if k2 not in scoresblast.keys(): scoresblast[k2] = {"score": 0, "name": "", "quality": 0.0}
                    score = blastScore(allalignments[k1][k2][0], allalignments[k1][k2][1])
                    if score > scoresblast[k2]["score"]:
                        scoresblast[k2]["score"] = score
                        scoresblast[k2]["name"] = k1
                        scoresblast[k2]["quality"] = blastQuality(allalignments[k1][k2][0], allalignments[k1][k2][1])

    return scoresblast

if __name__ == '__main__':
    clocal = False
    cglobal = False
    miniblast = False
    blastbonus = False
    processblast = True

    accession_codes = {
        # 7 known human coronaviruses
        "Human-SARS-CoV-2": "NC_045512",
        "Human-SARS": "NC_004718",
        "Human-MERS": "NC_019843",
        "Human-HCoV-OC43": "NC_006213",
        "Human-HCoV-229E": "NC_002645",
        "Human-HCoV-NL63": "NC_005831",
        "Human-HCoV-HKU1": "NC_006577",

        # Bat
        "Bat-CoV MOP1": "EU420138",
        "Bat-CoV HKU8": "NC_010438",
        "Bat-CoV HKU2": "NC_009988",
        "Bat-CoV HKU5": "NC_009020",
        "Bat-CoV RaTG13": "MN996532",
        "Bat-CoV-ENT": "NC_003045",

        # Other animals
        "Hedgehog-CoV 2012-174/GER/2012": "NC_039207",
        "Pangolin-CoV MP789": "MT121216",
        "Rabbit-CoV HKU14": "NC_017083",
        "Duck-CoV isolate DK/GD/27/2014": "NC_048214",
        "Feline infectious peritonitis virus": "NC_002306",  # cat
        "Giraffe-CoV US/OH3/2003": "EF424623",
        "Murine-CoV MHV/BHKR_lab/USA/icA59_L94P/2012": "KF268338",  # mouse
        "Equine-CoV Obihiro12-2": "LC061274",  # horse
    }

    extraction = {}
    for code in accession_codes:
        print(code, accession_codes[code])
        extraction[code] = extractSpikeProtein(accession_codes[code])[0]


    if clocal: computeLocal(extraction)
    if cglobal: computeGlobal(extraction)

    # allalignments = loadAlignments("alignments.txt")
    # sample = allalignments["Human-SARS-CoV-2"]["Human-SARS"]
    # global_alignment(sample[0], sample[1], scoring)
    # dendrogram(allalignments, "problem2.png")

    reference_genomes = [
        "Human-SARS",
        "Bat-CoV RaTG13",
        "Pangolin-CoV MP789",
    ]
    query = "Human-SARS-CoV-2"
    orf_candidates = {
        "ORF-1": (1, 11995, 13483),
        "ORF-2": (1, 26792, 27191),
        "ORF-3": (1, 23650, 25384),
        "ORF-4": (1, 9133, 13483),
        "ORF-5": (1, 25392, 26220),
    }

    if miniblast: miniBLAST(reference_genomes, query, orf_candidates)

    reference_genomes2 = [
        "Human-MERS",
        "Bat-CoV HKU5",
        "Hedgehog-CoV 2012-174/GER/2012",
    ]
    query2 = "Human-SARS-CoV-2"

    orf_candidates2 = {
        "ORF-1": (1, 11995, 13483),
        "ORF-2": (1, 26792, 27191),
        "ORF-3": (1, 23650, 25384),
        "ORF-4": (1, 9133, 13483),
        "ORF-5": (1, 25392, 26220),
    }

    if blastbonus: miniBLAST(reference_genomes2, query2, orf_candidates2, file="BONUS")

    if processblast:
        allalignments = loadAlignments("miniBlast.txt")
        rs = processBLASTResults(allalignments, orf_candidates)
        for r in rs: print(r, rs[r])
        dendrogram(allalignments, "miniblast.png", "!H")

        allalignments = loadAlignments("miniBlastBONUS.txt")
        rs = processBLASTResults(allalignments, orf_candidates2)
        for r in rs: print(r, rs[r])
        # dendrogram(allalignments, "miniblastBonus.png", "!H")