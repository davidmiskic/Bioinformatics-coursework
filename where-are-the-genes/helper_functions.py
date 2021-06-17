import pickle
from os import path
from typing import Tuple, Generator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import matplotlib.pyplot as plt

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


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[loc.start.position:loc.end.position]
        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, loc.start.position, loc.end.position))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((region.strand, loc.start.position, loc.end.position))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (loc.start.position, loc.end.position, region.strand)
                )
                print("\t", str(ex))

    return orfs


def isNotNested(start, end, others):
    for orf in others:
        if orf[0] <= start and orf[1] >= end: return False
    return True

def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

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

def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

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

def orf_to_seq(orf, record):
    if orf[0] == 1: return str(record.seq[orf[1]:orf[2]])
    else: str((record.seq[orf[1]:orf[2]]).reverse_complement())

def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """
    # TODO
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



def metrics(record, start_codons, stop_codons, interval = 9):
    orfs = set(extract_gt_orfs(record, start_codons, stop_codons))
    candidates = find_all_orfs(record.seq, start_codons, stop_codons)
    precs = []
    recs = []
    fs = []
    maxorf = int(max([abs(candidate[1] - candidate[2]) for candidate in candidates])) + 1
    for cutoff in range(0, maxorf, interval):
        candidatesCutoff = [candidate for candidate in candidates if abs(candidate[1] - candidate[2]) >= cutoff]
        #print(len(candidatesCutoff))
        candidatesCutoff = set(candidatesCutoff)
        TP = candidatesCutoff.intersection(orfs)
        FP = candidatesCutoff.difference(orfs)
        FN = (set(candidates).difference(candidatesCutoff)).intersection(orfs)
        precision = len(TP) / (len(TP)+len(FP))
        recall = len(TP) / (len(TP)+len(FN))
        #print(len(TP), len(FP), len(FN))
        #print(precision, recall)
        precs.append(precision)
        recs.append(recall)
        if (precision + recall) != 0: fs.append(2 * (precision * recall) / (precision + recall))
        else: fs.append(0)

    together = [precs[i]+recs[i]+fs[i] for i in range(0, len(fs))]
    # which interval does index represent
    optimalCutoff = together.index(max(together))*interval
    plt.plot([x for x in range(0, maxorf, interval)], recs, color="blue", label="Recall")
    plt.plot([x for x in range(0, maxorf, interval)], precs, color="red", label="Precision")
    plt.plot([x for x in range(0, maxorf, interval)], fs, color="green", label="F1")
    # plt.axis([0, 1, 0, 1])
    plt.xlabel('cutoff[#bases]')
    plt.legend()
    plt.savefig('plot.png', dpi=250)
    plt.show()
    return optimalCutoff


#record = load("NC_006058")
Entrez.email = "<your email>@fri.uni-lj.si"
#start_codons = ["ATG"]
#stop_codons = ["TGA"]
#metrics(record, start_codons, stop_codons)

# b = set(extract_gt_orfs(record, start_codons, stop_codons))
# b1 = [orf for orf in b if orf[0] == -1]
# b2 = [orf for orf in b if orf[0] == 1]
# a = find_all_orfs(record.seq, start_codons, stop_codons)
# a1 = [orf for orf in a if orf[0] == -1]
# a2 = [orf for orf in a if orf[0] == 1]
# c1 = set(a1).intersection(set(b1))
# c2 = set(a2).intersection(set(b2))
# print("-", len(a1), len(b1), len(c1))
# print("+", len(a2), len(b2), len(c2))



# print("{:,} bases".format(len(record.seq)))
# print("Finding ORFs using start/stop codons...")
# orf_candidates = find_all_orfs(record.seq, start_codons, stop_codons)
# print(f"{len(orf_candidates)} ORF candidates found")