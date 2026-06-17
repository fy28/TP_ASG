import streamlit as st
import io
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter
from Bio.SeqRecord import SeqRecord

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


uploaded_file = st.file_uploader(
    "Choisissez un fichier FASTQ",
    type=["fastq"]
)

if uploaded_file is not None:

    # binaire => texte
    text_file = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

    # SeqIO.parse lit le fichier FASTQ et renvoie les reads un par un
    # resultat => list pour compter facilement le nombre total de reads
    reads = list(SeqIO.parse(text_file, "fastq"))

    st.success("Fichier chargé avec succès !")

    st.write(f"Nombre total de reads : **{len(reads)}**")

    st.subheader("Aperçu des premiers reads")

        # parcour des 5 premiers reads de la liste, si < 5, arret
    for read in reads[:5]:

        st.write(f"ID : {read.id}")

        # séquence => texte
        sequence = str(read.seq)

        st.code(sequence)

        kmers = generate_kmers(sequence, k)

        # aff
        st.write(f"Nombre de k-mers : {len(kmers)}")

        st.write(kmers)

        # Séparation visuelle entre deux reads
        st.divider()

    st.subheader("Histogramme des fréquences des k-mers")

    kmer_counts = count_kmers(reads, k)

    fig, ax = plt.subplots()

    # Construction de l'histogramme des fréquences
    ax.hist(
        kmer_counts.values(),
        bins=range(1, max(kmer_counts.values()) + 2),
    )

    #
    ax.set_xlabel("Fréquence")
    ax.set_ylabel("Nombre de k-mers")
    ax.set_title("Distribution des fréquences des k-mers")

    # Aff
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