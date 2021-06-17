import sys
sys.setrecursionlimit(1500)

def kmers(seq, k, stride=1):
    """Generate a list of k-mers from a given string.

    Parameters
    ----------
    seq: str
        A string which will be decomposed into k-mers.
    k: int
    stride: int

    Returns
    -------
    List[str]

    Examples
    --------
    >>> kmers("mesenchyme", k=3, stride=1)
    ["mesench", "esenchy", "senchym", "enchyme"]

    >>> kmers("mesenchyme", k=3, stride=2)
    ["mesench", "senchym"]

    """
    ind = 0
    kmersList = []
    while (ind <= len(seq) - k):
        kmersList.append(seq[ind:ind + k])
        ind += stride
    return kmersList


# print(kmers("mesenchyme", 7, 2))
# print(kmers("ABC", 2, 1))

def assemble_genome(seqs, k=None, stride=None):
    """Perform genome assembly using the Eulerian path algorithm.

    The overall algorithm should follow the following general structure:
    1. For an input list of sequences, e.g. kmers, construct a DeBruijn graph as
       seen in the lectures.
    2. Find all possible Euerlian paths through the graph, i.e. all possible paths
       which visit each edge exactly once. Your paths should all start from a
       source node with in-degree zero. In case no such node exists, you may use
       the first sequence in the list of input sequences as your starting point.
    3. Decode your obtained paths into sequences, and return a list (or set) of
       unique genome assemblies as strings.

    Parameters
    ----------
    seqs: List[str]
        The list of strings. In this homework, this will always be a list of k-mers.

    k: int, optional
        This parameter is optional, and you may ignore it if you do not need it.
        But, this function will be called with the same `k` as `kmers`, which will
        be used to produce `seqs`. Therefore, make sure you do not remove this
        parameter.
    stride: int, optional
        This parameter is optional, and you may ignore it if you do not need it.
        But, this function will be called with the same `stride` as `kmers`, which will
        be used to produce `seqs`. Therefore, make sure you do not remove this
        parameter.

    Returns
    -------
    Set[str]
        A set (or list if you really want) of unique assemblies for the given `seqs`.

    """
    nodes = set()
    degrees = []
    edges = []  # [from, to, label]
    # make de Bruijn graph from k-mers
    for seq in seqs:
        oneLess = kmers(seq, k - 1, 1)
        i = 0
        for node in oneLess: nodes.add(node)
        while (i < len(oneLess) - 1):
            edges.append([oneLess[i], oneLess[i + 1], seq])
            i += 2
    # check if graph is Eulerian

    unbalanced = []
    inDegZero = []
    for node in nodes:
        degNodeOut = 0
        degNodeIn = 0
        for edge in edges:
            if edge[0] == node: degNodeOut += 1
            if edge[1] == node: degNodeIn += 1
        degrees.append([node, degNodeIn, degNodeOut])
        if degNodeIn == 0: inDegZero.append([node, degNodeIn, degNodeOut])
        if degNodeOut != degNodeIn: unbalanced.append([node, degNodeIn, degNodeOut])
    if len(inDegZero) > 0:
        #print("A-DEGZERO", inDegZero)
        temp = inDegZero[0]
        for i in range(0, len(edges)):
            if edges[i][0] == temp[0]:
                newE = edges.pop(i)
                break
        findCycle3(edges, newE, len(edges) + 1, [newE[0]])
        edges.append(newE)
    else:
        findCycle3(edges[1:], edges[0], len(edges), [edges[0][0]])

    global paths
    founds = []
    output = set()
    for path in paths:
        if len(path) == len(edges) + 1 and path not in founds:
            founds.append(path)
            merge = [path[0]]
            for p in path[1:]: merge.append(p[-1])
            merge = "".join(merge)
            output.add(merge)
    paths = []
    return output

# store found paths found in recursive calls
paths = []

def findCycle3(edges, startEdge, nedges, foundPath):
    # eulerjev cikel: obiščeš vse povezave
    if len(edges) == 0 and len(foundPath) == nedges + 1:
        global paths
        paths.append(foundPath)
        return foundPath
    if len(foundPath) < 2:
        path = foundPath + [startEdge[1]]
    else:
        path = foundPath
    start = startEdge[1]
    for i in range(0, len(edges)):
        edge = edges[i]
        if edge[0] == start:
            # če je pot od trenutnega do naslednjega, pojdi po tej poti in odstrani iz možnih
            # trenutni je konec te poti
            newL = edges.copy()
            newStart = newL.pop(i)
            npath = path.copy()
            npath.append(edge[1])
            if len(foundPath) < nedges + 2:
                findCycle3(newL, newStart, nedges, npath)
            else:
                return

# seqs_test = kmers("AAABBBA", 3, 1)
# assemble_genome(seqs_test, 3, 1)
#
# a = seqs_test = kmers("to_every_thing_turn_turn_turn_there_is_a_season", 4, 1)
# assemble_genome(seqs_test, 4, 1)
#
# seqs_test = kmers("DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDGVVFQNRESVLPTQSYGYYHEYTVITPGARTRGTRRIITGEATQEDYYTGDHYATFSLIDQTC", 3, 1)
# assemble_genome(seqs_test, 3, 1)

spikeProtein = "MFLLTTKRTMFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVY" \
               "FASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLRE" \
               "FVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDC" \
               "ALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLC" \
               "FTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFP" \
               "LQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCS" \
               "FGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQS" \
               "IIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKT" \
               "PPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAG" \
               "AALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEV" \
               "QIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPR" \
               "EGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAK" \
               "NLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

#
#
# for i in range(10, 2, -1):
#     minimum_fragment_length = i
#     a = assemble_genome(kmers(spikeProtein, minimum_fragment_length, 1), minimum_fragment_length, 1)
#     print(i, len(a))