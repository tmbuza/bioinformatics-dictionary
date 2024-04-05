import csv
import json

# Existing bioinformatics dictionary
bioinformatics_dictionary = {
    "NGS": "Next-Generation Sequencing",
    "SNPs": "Single Nucleotide Polymorphisms",
    "Alignment": "Process of aligning sequencing reads to a reference genome",
    "Variant Calling": "Identification of genetic variations such as SNPs and indels",
    "Differential Expression Analysis": "Comparison of gene expression levels between conditions",
    "Variant Annotation": "Adding functional annotations to identified genetic variants",
    # Add more terms as needed
}

# New terms to add
new_terms = {
    "FASTQ": "File format containing raw sequencing reads and their quality scores",
    "Variant Filtering": "Process of filtering genetic variants based on quality metrics",
    # Add more new terms as needed
}

# Update the dictionary with new terms
for term, definition in new_terms.items():
    if term not in bioinformatics_dictionary:
        bioinformatics_dictionary[term] = definition
    else:
        # Concatenate the definitions if the term already exists
        bioinformatics_dictionary[term] += f", {definition}"

# Convert dictionary to list of lists (rows) for CSV
rows = [["Term", "Definition"]]
for term, definition in bioinformatics_dictionary.items():
    rows.append([term, definition])

# Write data to CSV file
with open('dictionary.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(rows)

# Existing bioinformatics dictionary
bioinformatics_dictionary = {
    "NGS": "Next-Generation Sequencing",
    "SNPs": "Single Nucleotide Polymorphisms",
    "Alignment": "Process of aligning sequencing reads to a reference genome",
    "Variant Calling": "Identification of genetic variations such as SNPs and indels",
    "Differential Expression Analysis": "Comparison of gene expression levels between conditions",
    "Variant Annotation": "Adding functional annotations to identified genetic variants",
    # Add more terms as needed
}

# New terms to add
new_terms = {
    "FASTQ": "File format containing raw sequencing reads and their quality scores",
    "Variant Filtering": "Process of filtering genetic variants based on quality metrics",
    # Add more new terms as needed
}

# Update the dictionary with new terms
for term, definition in new_terms.items():
    if term not in bioinformatics_dictionary:
        bioinformatics_dictionary[term] = definition
    else:
        # Concatenate the definitions if the term already exists
        bioinformatics_dictionary[term] += f", {definition}"

# Convert dictionary to HTML
html_content = "<html><head><title>Bioinformatics Dictionary</title></head><body>"
html_content += "<h1>Bioinformatics Dictionary</h1>"
html_content += "<table border='2'><tr><th>Term</th><th>Definition</th></tr>"
for term, definition in bioinformatics_dictionary.items():
    html_content += f"<tr><td>{term}</td><td>{definition}</td></tr>"
html_content += "</table></body></html>"

# Write HTML content to file
with open('dictionary.html', 'w') as html_file:
    html_file.write(html_content)
