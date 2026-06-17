import streamlit as st
import io
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter
from bitarray import bitarray
import mmh3

st.title("Projet ASG-2026")
st.header("Lot 1 : Lecture de fichiers FASTQ")

# Choix de la taille des kmers
k = st.slider(
    "Choisissez la taille des k-mers",
    min_value=2,
    max_value=15,
    value=3
)


def generate_kmers(sequence, k):
    """
    Génère tous les k-mers d'une séquence.

    Exemple :
    sequence = "ATCGAT"
    k = 3

    Retourne :
    ["ATC", "TCG", "CGA", "GAT"]
    """

    kmers = []

    # Parcourir la séquence en avançant d'un nucléotide à la fois
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)

    return kmers


def count_kmers(reads, k):
    """
    Compte la fréquence de tous les k-mers.
    """

    all_kmers = []

    for read in reads:
        sequence = str(read.seq)
        kmers = generate_kmers(sequence, k)
        all_kmers.extend(kmers)

    return Counter(all_kmers)


# LCS = Longest Common Subsequence
def longest_common_subsequence(seq1, seq2):

    m = len(seq1)
    n = len(seq2)

    # Matrice de programmation dynamique
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Remplissage de la matrice
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1

            else:
                dp[i][j] = max(
                    dp[i - 1][j],
                    dp[i][j - 1]
                )

    # Reconstruction de la sous-séquence
    lcs = []

    i = m
    j = n

    while i > 0 and j > 0:

        if seq1[i - 1] == seq2[j - 1]:
            lcs.append(seq1[i - 1])
            i -= 1
            j -= 1

        elif dp[i - 1][j] > dp[i][j - 1]:
            i -= 1

        else:
            j -= 1

    lcs.reverse()

    #return "".join(lcs), dp[m][n]
    #
    lcs_string = "".join(lcs)

    start_seq1 = seq1.find(lcs_string)
    start_seq2 = seq2.find(lcs_string)

    return (
        lcs_string,
        dp[m][n],
        start_seq1,
        start_seq2
    )

# Bloom Filter
class BloomFilter:

    def __init__(self, size=1000, hash_count=3):

        self.size = size
        self.hash_count = hash_count

        # tableau de bits
        self.bit_array = bitarray(size)
        self.bit_array.setall(0)

    # insertion d'un k-mer
    def add(self, item):

        for i in range(self.hash_count):

            index = mmh3.hash(
                item,
                i
            ) % self.size

            self.bit_array[index] = 1

    # test de présence
    def contains(self, item):

        for i in range(self.hash_count):

            index = mmh3.hash(
                item,
                i
            ) % self.size

            if not self.bit_array[index]:
                return False

        return True


def build_bloom_filter(reads, k):

    bloom = BloomFilter()

    for read in reads:

        sequence = str(read.seq)

        kmers = generate_kmers(
            sequence,
            k
        )

        for kmer in kmers:
            bloom.add(kmer)

    return bloom

#k-mers neighbor
def get_successors(kmer, bloom):
    """
    Génère les 4 extensions possibles en 3'.
    """

    bases = ["A", "C", "G", "T"]
    successors = []

    for base in bases:
        next_kmer = kmer[1:] + base

        if bloom.contains(next_kmer):
            successors.append(next_kmer)

    return successors


uploaded_file = st.file_uploader(
    "Choisissez un fichier FASTQ",
    type=["fastq"]
)

if uploaded_file is not None:

    # binaire => texte
    text_file = io.StringIO(
        uploaded_file.getvalue().decode("utf-8")
    )

    # SeqIO.parse lit le FASTQ et renvoie les reads
    reads = list(SeqIO.parse(text_file, "fastq"))

    st.success("Fichier chargé avec succès !")

    st.write(f"Nombre total de reads : **{len(reads)}**")

    st.subheader("Aperçu des premiers reads")

    # parcour des 5 premiers reads de la liste
    for read in reads[:5]:

        st.write(f"ID : {read.id}")

        # séquence => texte
        sequence = str(read.seq)

        st.code(sequence)

        kmers = generate_kmers(sequence, k)

        st.write(f"Nombre de k-mers : {len(kmers)}")

        st.write(kmers)

        # séparation visuelle
        st.divider()

    st.subheader("Histogramme des fréquences des k-mers")

    # Comptage global des k-mers
    kmer_counts = count_kmers(reads, k)

    fig, ax = plt.subplots()

    # Construction de l'histogramme
    ax.hist(
        kmer_counts.values(),
        bins=range(
            1,
            max(kmer_counts.values()) + 2
        )
    )

    ax.set_xlabel("Fréquence")
    ax.set_ylabel("Nombre de k-mers")
    ax.set_title("Distribution des fréquences des k-mers")

    st.pyplot(fig)

    # ===== FASTQ -> FASTA =====

    st.subheader("Conversion FASTQ → FASTA")

    fasta_content = ""

    for read in reads:
        fasta_content += f">{read.id}\n{read.seq}\n"

    st.download_button(
        label="Télécharger le fichier FASTA",
        data=fasta_content,
        file_name="converted.fasta",
        mime="text/plain"
    )

# ======================
# LOT 2
# ======================

st.header("Lot 2 : Alignement")

seq1 = st.text_input(
    "Read 1",
    value="ATCGAT"
)

seq2 = st.text_input(
    "Read 2",
    value="TCGATA"
)

if seq1 and seq2:

    # Calcul du LCS
    #old one
    # lcs, score = longest_common_subsequence(
    #     seq1,
    #     seq2
    # )
    lcs, score, pos1, pos2 = longest_common_subsequence(
    seq1,
    seq2
    )

    st.write(f"Longueur Read 1 : {len(seq1)}")
    st.write(f"Longueur Read 2 : {len(seq2)}")

    st.write(f"Score : {score}")
    #barre visuelle de similarite mais pas oblig
    st.progress(score / max(len(seq1), len(seq2)))
    st.write(f"Position dans Read 1 : {pos1}")

    st.write(f"Position dans Read 2 : {pos2}")

    st.write(
        f"Sous-séquence commune : {lcs}"
    )

# ======================
# LOT 3
# ======================

st.header("Lot 3 : Bloom Filter")

if uploaded_file is not None:

    # Construction du filtre
    bloom = build_bloom_filter(
        reads,
        k
    )

    kmer_test = st.text_input(
        "Tester un k-mer"
    )

    if kmer_test:

        if bloom.contains(kmer_test):

            st.success(
                "Probablement présent"
            )

        else:

            st.error(
                "Certainement absent"
            )

    st.subheader("Voisins 3' du k-mer")

    if kmer_test:

        successors = get_successors(
            kmer_test,
            bloom
        )

        st.write("Voisins trouvés :")
        st.write(successors)

        if len(successors) == 0:
            st.warning("Aucun voisin trouvé")

        elif len(successors) == 1:
            st.success("Progression unique possible")

        else:
            st.warning("Bifurcation détectée")