import requests
import xml.etree.ElementTree as ET
import sys

# List of UniProt accession IDs for the proteins you want to retrieve
proteins = open(sys.argv[1])

while True:
    protein = proteins.readline().strip()
    if not protein:
        break

    # Define the URL of the XML data for each protein
    url = f'https://rest.uniprot.org/uniprotkb/{protein}.xml'

    # Send an HTTP GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (HTTP status code 200)
    if response.status_code == 200:
        # Parse the XML data
        xml_data = response.content

        # Parse the XML content
        root = ET.fromstring(xml_data)

        # Find and retrieve the <fullName> field
        full_name = None
        for entry in root.findall(".//{http://uniprot.org/uniprot}fullName"):
            full_name = entry.text
            break  # Stop after the first name is found

        # If no name is retrieved, set it to "Predicted Protein"
        if not full_name:
            full_name = "Predicted Protein"

        # Find and retrieve the <name> field within <gene> section
        gene_name = None
        for gene_entry in root.findall(".//{http://uniprot.org/uniprot}gene/{http://uniprot.org/uniprot}name[@type='primary']"):
            gene_name = gene_entry.text
            break  # Stop after the first gene name is found

        # If no gene name is retrieved, set it to "Gene Name Not Found"
        if not gene_name:
            gene_name = "Gene Name Not Found"

        # Print the full name and gene name
        print(f'{protein}\tprotein:{protein} gene_symbol:{gene_name} description:{full_name}') 

    else:
        print(f"Failed to retrieve data for {protein}. HTTP status code:", response.status_code)

