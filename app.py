import streamlit as st
import io
from Bio import SeqIO

st.title("Projet ASG-2026")
st.header("Lot 1 : Lecture de fichiers FASTQ")

# Choix de la taille des k-mers par l'utilisateur
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

        # Conversion séquence => texte
        sequence = str(read.seq)

        st.code(sequence)

        # Génération des k-mers à partir de la séquence
        kmers = generate_kmers(sequence, k)

        # aff
        st.write(f"Nombre de k-mers : {len(kmers)}")

        # Aff liste des k-mers
        st.write(kmers)

        # Séparation visuelle entre deux reads
        st.divider()