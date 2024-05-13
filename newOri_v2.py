import io
import streamlit as st
import requests
from Bio import SeqIO, SeqUtils
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import time
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# KD dictionary
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
      "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
      "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Function to fetch protein data from UniProt
def fetch_protein_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    record = SeqIO.read(io.StringIO('\n'.join(response.text.splitlines())), "fasta")

    try:
        molecular_weight = SeqUtils.molecular_weight(record.seq)
    except ValueError:
        molecular_weight = SeqUtils.molecular_weight(record.seq, seq_type="protein")

    return {
        "sequence": str(record.seq),
        "length": len(record.seq),
        "molecular_weight": molecular_weight
    }


# Function to fetch protein-protein interaction network from STRING DB
# Function to fetch protein-protein interaction network from STRING DB
def fetch_ppi_network(uniprot_id):
    url = f"https://string-db.org/api/json/interaction_partners?identifiers={uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for 4xx and 5xx status codes
        data = response.json()

        # Check if data is a list of dictionaries
        if isinstance(data, list):
            return data
        else:
            raise ValueError("Unexpected data structure")

    except Exception as e:
        st.error("An unexpected error occurred while fetching protein-protein interaction network:", e)
        return None

# Function to visualize protein-protein interaction network
def visualize_ppi_network(data):
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes and edges to the graph
    for interaction in data:
        protein_A = interaction["preferredName_A"]
        protein_B = interaction["preferredName_B"]
        G.add_edge(protein_A, protein_B)

    # Draw the graph
    fig, ax = plt.subplots()
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, edge_color='black', linewidths=1, font_size=10)
    plt.title("Protein-Protein Interaction Network")
    st.pyplot(fig)

# Function to perform sequence alignment
def perform_sequence_alignment(protein_sequence):
    # URL for the Clustal Omega REST API
    url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"

    # Parameters for the API request
    params = {
        "sequence": protein_sequence,
        # "email": "your-email@example.com",  # Replace with your email
    }

    # Send a POST request to the API
    response = requests.post(url, data=params)
    if not response.ok:
        print(f"Failed to start alignment. HTTP status code: {response.status_code}")
        print(f"Response text: {response.text}")
        return None

    # Get the job ID from the response
    job_id = response.text

    # URL to get the status of the job
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"

    # Wait for the job to finish
    while True:
        response = requests.get(status_url)
        if response.text == "FINISHED":
            break
        time.sleep(1)

    # URL to get the results of the job
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal_num"

    # Get the results
    response = requests.get(result_url)
    if not response.ok:
        print(f"Failed to get alignment results. HTTP status code: {response.status_code}")
        print(f"Response text: {response.text}")
        return None

    # Return the aligned sequence
    return response.text

def main():
    st.title("Protein Data Analysis")

    option = st.sidebar.selectbox(
        'Choose Input Type',
        ('UniProt ID', 'Protein Sequence'))

    if option == 'UniProt ID':
        uniprot_id = st.sidebar.text_input('Enter UniProt ID')

        if st.sidebar.button('Fetch Data'):
            protein_data = fetch_protein_data(uniprot_id)
            if protein_data:
                st.write("### Protein Characteristics")
                st.write(f"Length: {protein_data['length']}")
                st.write(f"Molecular Weight: {protein_data['molecular_weight']}")

                st.write("### Protein-Protein Interaction Network")
                ppi_network = fetch_ppi_network(uniprot_id)
                if ppi_network:
                    visualize_ppi_network(ppi_network)
                else:
                    st.write("Failed to fetch PPI network.")

    elif option == 'Protein Sequence':
        protein_sequence = st.sidebar.text_area('Enter Protein Sequence')

        if st.sidebar.button('Analyze Sequence'):
            st.write("### Protein Characteristics")
            length = len(protein_sequence)
            molecular_weight = SeqUtils.molecular_weight(protein_sequence, seq_type="protein")
            st.write(f"Length: {length}")
            st.write(f"Molecular Weight: {molecular_weight}")

            st.write("### Sequence Alignment")
            alignment_output = perform_sequence_alignment(protein_sequence)
            st.write(alignment_output)

            st.write("### Additional Merits")
            protein_analysis = ProteinAnalysis(protein_sequence)
            isoelectric_point = protein_analysis.isoelectric_point()
            st.write(f"Isoelectric Point: {isoelectric_point}")

            hydrophobicity = protein_analysis.protein_scale(param_dict=KD, window=9, edge=1.0)[0]
            st.write(f"Hydrophobicity: {hydrophobicity}")

main()
