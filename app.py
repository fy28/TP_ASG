import streamlit as st

# Titre de la page
st.title("🧬 Projet ASG-2026")

# Petit texte de bienvenue
st.write("Bienvenue dans notre application d'assemblage de séquences !")

# Champ de saisie
nom = st.text_input("Quel est ton prénom ?")

# Afficher un message si quelque chose est saisi
if nom:
    st.success(f"Bonjour {nom} ! 👋")