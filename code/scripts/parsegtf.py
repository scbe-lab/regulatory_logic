import argparse
import re

# Define regular expressions for matching the fields in column 9
gene_id_re = re.compile('gene_id=([^;]+)')
transcript_id_re = re.compile('ID=([^;]+)')
parent_re = re.compile('Parent=([^;]+)')

# Parse the command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='path to the input GTF file')
parser.add_argument('output_file', help='path to the output GTF file')
args = parser.parse_args()

# Open the input and output files
with open(args.input_file, 'r') as f_in, open(args.output_file, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            continue
        # Split the line into fields
        fields = line.strip().split('\t')
        # Parse the fields in column 9
        col9_fields = fields[8].split(';')
        col9_dict = {}
        for field in col9_fields:
            key, value = field.strip().split('=')
            col9_dict[key] = value.strip('"')
        # If the row corresponds to a gene
        if fields[2] == 'gene':
            # Get the gene ID and add it to the column 9 dictionary
            gene_id_match = gene_id_re.search(fields[8])
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                col9_dict['gene_name'] = gene_id
        # If the row corresponds to a transcript
        elif fields[2] == 'transcript':
            # Get the transcript ID and gene ID and add them to the column 9 dictionary
            transcript_id_match = transcript_id_re.search(fields[8])
            gene_id_match = gene_id_re.search(fields[8])
            if transcript_id_match and gene_id_match:
                transcript_id = transcript_id_match.group(1)
                gene_id = gene_id_match.group(1)
                col9_dict['gene_name'] = gene_id
                col9_dict['transcript_name'] = transcript_id
        # If the row corresponds to an exon, CDS, 5' UTR, or 3' UTR
        else:
            # Get the parent ID and add it to the column 9 dictionary
            parent_match = parent_re.search(fields[8])
            if parent_match:
                parent_id = parent_match.group(1)
                # Remove any last termination preceded by a `.`
                parent_id = parent_id.rsplit('.', 1)[0]
                col9_dict['gene_name'] = parent_id
                col9_dict['transcript_name'] = parent_match.group(1)
        # Reconstruct the column 9 string
        col9_fields = []
        for key, value in col9_dict.items():
            col9_fields.append('{} "{}"'.format(key, value))
        col9_string = '; '.join(col9_fields)
        # Write the modified line to the output file
        fields[8] = col9_string
        f_out.write('\t'.join(fields) + '\n')