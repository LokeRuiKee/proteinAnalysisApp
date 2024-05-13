# Import necessary libraries
import streamlit as st
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import requests
import pandas as pd
from bs4 import BeautifulSoup
import networkx as nx
import matplotlib.pyplot as plt

# Define function to retrieve protein data from Uniprot
def fetch_protein_data(uniprot_id):
    # Make API request to Uniprot
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        return None
    
# Define function to retrieve protein-protein interaction network from STRING DB

# convert uniprot_id to protein name
def uniprot_id_to_protein_name(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the XML response to get the protein name
        soup = BeautifulSoup(response.content, 'xml')
        protein_name = soup.find('fullName').text
        return protein_name
    else:
        return None

def create_graph_from_ppi_data(ppi_data):
    # Create a graph from the PPI data
    ppi_graph = nx.from_pandas_edgelist(ppi_data, "preferredName_A", "preferredName_B")
    return ppi_graph

def visualize_graph(ppi_graph):
    # Visualize the PPI graph
    pos = nx.spring_layout(ppi_graph)
    nx.draw(ppi_graph, pos, with_labels=True)
    plt.show()

def retrieve_ppi(uniprot_id):
    protein_name = uniprot_id_to_protein_name(uniprot_id)
    
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": protein_name,
        "species": 9606  # Homo sapiens
    }

    response = requests.get(string_url, params=params)
    data = response.json()
    
    network_df = pd.json_normalize(data)
    
    return network_df

# PPI network
def visualize_ppi(protein_name):
    # Retrieve the PPI data for the given protein
    ppi_data = retrieve_ppi(protein_name)
    
    # Create a graph from the PPI data
    ppi_graph = create_graph_from_ppi_data(ppi_data)
    
    # Visualize the PPI graph
    visualize_graph(ppi_graph)

# Define function to perform sequence alignment
def perform_sequence_alignment(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments

# Main function to run the Streamlit web app
def main():
    st.title("Protein Data Analysis App")
    
    # Sidebar for user input
    option = st.sidebar.selectbox(
        'Select Input Type',
        ('Uniprot ID', 'Protein Sequence')
    )

    if option == 'Uniprot ID':
        uniprot_id = st.sidebar.text_input('Enter Uniprot ID:')
        if uniprot_id:
            protein_data = fetch_protein_data(uniprot_id)
            if protein_data:
                # Display protein characteristics
                st.header("Protein Characteristics")
                st.write(protein_data)
                
                # Display protein-protein interaction network
                st.header("Protein-Protein Interaction Network")
                visualize_ppi(uniprot_id)
                if visualize_ppi:
                    st.write(visualize_ppi)
                else:
                    st.write("No interaction network found.")
            else:
                st.write("Invalid Uniprot ID.")
    
    elif option == 'Protein Sequence':
        protein_sequence = st.sidebar.text_area('Enter Protein Sequence:')
        if protein_sequence:
            # Perform sequence alignment
            st.header("Sequence Alignment")
            seq2 = "Example sequence"  # You need to provide an example sequence for alignment
            alignments = perform_sequence_alignment(protein_sequence, seq2)
            for alignment in alignments:
                st.text(format_alignment(*alignment))

main()