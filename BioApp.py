# Import necessary libraries
import streamlit as st
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import requests
import pandas as pd

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
def fetch_ppi_network(uniprot_id):
    url = f"https://string-db.org/api/tsv/network?identifier={uniprot_id}"
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the TSV response to get the PPI network data
        ppi_network_data = pd.read_csv(pd.compat.StringIO(response.text), sep="\t")
        return ppi_network_data
    else:
        return None

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
                ppi_network = fetch_ppi_network(uniprot_id)
                if ppi_network is not None:
                    st.write(ppi_network)
                else:
                    st.write("No PPI network data available.")
# Run the app
main()