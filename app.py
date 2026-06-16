import streamlit as st
import io
from Bio import SeqIO

st.title(" Projet ASG-2026")
st.header("Lot 1 : Lecture de fichiers FASTQ")

uploaded_file = st.file_uploader(
    "Choisissez un fichier FASTQ",
    type=["fastq"]
)

if uploaded_file is not None:
    # binaire => texte 
    text_file = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

    reads = list(SeqIO.parse(text_file, "fastq"))

    st.success(f"Fichier chargé avec succès !")

    st.write(f"Nombre total de reads : **{len(reads)}**")

    st.subheader("Aperçu des premiers reads")

    #parcour des 5 premiers reads de la liste, si < 5, arret 
    for read in reads[:5]:
        st.write(f"ID : {read.id}")
        st.code(str(read.seq))

#SeqIO.parse lit le fichier FASTQ et renvoie les reads un par un
#resultat => list pour compter facilement le nombre total de reads