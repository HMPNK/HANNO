# Author: C. Dassio
#
# This program takes a bedDB (created by HANNO) or BED12 file and produces a GTF or GFF3 file (with optional information chosen by the user) (default GTF)
# v1: supports GFF format for bedDB and BED12 input
# v2: added optional second input (genomic file) to give sequence for every entry
# v3: parsing the genome is now done with SeqIO from Biopython
# v4: parallel processing/generating of the GTF/GFF lines
# v5: biopython and parallel processing were discarded (took longer), some formating changes
# v6: the gene field (column 25) now has unique names, switched the transcript/mRNA sequence to concatenated exon and complete coding sequence, for sequences on the - strand, now gives the reverse complement
# v7: fixed formatting issues, supports GFF format for all BED inputs
# v8: removed BED3 functions, removed unique names for gene field, fixed formatting for sequence features in GFF, added sequence output for BED12
# v9: added option to create fasta files for mRNA and CDS sequences, fixed output-file-directory to be the same as the input-file directory if -o is given

import pandas as pd
import argparse
import sys
import os

# Helper functon to parse the genome fasta file
def parse_fasta(file_path):
    seq_dict = {}
    with open(file_path, 'r') as f:
        chrom = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if chrom:
                    seq_dict[chrom] = ''.join(sequence)
                # Start new sequence
                chrom = line[1:].split(' ')[0]  # Extract chromosome name from header
                sequence = []
            else:
                sequence.append(line)
        # Save the last sequence
        if chrom:
            seq_dict[chrom] = ''.join(sequence)
    return seq_dict

def load_genomic_sequences(genomic_file):
    return parse_fasta(genomic_file)

# Helper function to get the exact sequence from the genome data
def extract_sequence(seq_dict, chrom, start, end, strand):
    # Extract sequence snippet from the genomic sequences dictionary
    if chrom in seq_dict:
        sequence = seq_dict[chrom][start:end]  # Get the sequence slice
        if strand == '-':
            return reverse_complement(sequence)  # Get reverse complement if on negative strand
        return sequence  # Return sequence as is for positive strand
    else:
        return ""

# Helper function to give correct sequence depending on the strand
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}  # Dictionary for complement
    # Replace each base with its complement and reverse the sequence
    return ''.join(complement[base] for base in reversed(seq))

# Helper function to ensure correct gene entries
def collect_gene_positions(bedDB): 
    gene_positions = {}
    
    for _, row in bedDB.iterrows():
        mRNA_start = int(row.iloc[1])
        mRNA_end = int(row.iloc[2])
        gene_id = row.iloc[3].rsplit('.', 1)[0]
        
        if gene_id not in gene_positions:
            gene_positions[gene_id] = [mRNA_start, mRNA_end]
        else:
            current_start, current_end = gene_positions[gene_id]
            gene_positions[gene_id] = [min(current_start, mRNA_start), max(current_end, mRNA_end)]
    
    return gene_positions

# Helper function to format the sequences in the fasta files to have 80 characters per line
def write_fasta_with_formatting(file_path, header, sequence): 
    with open(file_path, 'a') as f:
        f.write(f"{header}\n")
        
        # Write the sequence in lines of 80 characters
        for i in range(0, len(sequence), 80):
            f.write(f"{sequence[i:i+80]}\n")

# Helper function to ensure fasta files are newly created/overwritten everytime the script is used and correct file path is used
def initialize_fasta_files(output_dir, output_file):
    # Get output base name without extension
    output_base_name = os.path.splitext(os.path.basename(output_file))[0]

    # Define FASTA file paths using the output file as a prefix
    mrna_fasta_path = os.path.join(output_dir, f"{output_base_name}_mRNA_sequences.fasta")
    cds_fasta_path = os.path.join(output_dir, f"{output_base_name}_CDS_sequences.fasta")

    # Overwrite the files by opening them in 'w' mode
    with open(mrna_fasta_path, 'w') as f:
        f.write("")  # Clear the file

    with open(cds_fasta_path, 'w') as f:
        f.write("")  # Clear the file

    # Return paths to be used later for appending
    return mrna_fasta_path, cds_fasta_path


# Function to convert one bedDB entry to GTF format
def bedDB_to_gtf(row, field_map, column_names, gene_positions, gene_ids_seen, seq_dict=None, include=None, fasta_enabled=False, input_file_path=None, mrna_fasta_path=None, cds_fasta_path=None):
    gtf_lines = []

    chrom = row.iloc[0]               # row.iloc[] was used to avoit problems with future pandas versions
    mRNA_start = int(row.iloc[1])
    mRNA_end = int(row.iloc[2])
    strand = row.iloc[5]
    cds_start = int(row.iloc[6] + 1)  # +1 due to GTF format being 1-indexed, unlike 0-indexed bed format
    cds_end = int(row.iloc[7])
    block_count = int(row.iloc[9])
    block_len = [int(x) for x in row.iloc[10].split(',') if x]
    block_starts = [int(x) for x in row.iloc[11].split(',') if x]

    # Check if exon sizes and starts match
    if len(block_len) != len(block_starts):
        raise ValueError("Problem with number of exons in bedDB")

    gene_id = row.iloc[3].rsplit('.', 1)[0]
    transcript_id = row.iloc[3]
    db_xref = row.iloc[3]
    attributes2 = ""
    attributes3 = ""
    attributes4 = ""

    # Handle --add_fields and --add_field_names
    if not field_map:
        attributes2 += f' gene "{row.iloc[25]}"; product "{row.iloc[14]}"; '
    else:
        added_fields = set()
        for index, name in field_map.items():
            if index < len(row):
                value = row.iloc[index]
                if name not in added_fields:
                    attributes3 += f'{name} "{value}"; '
                    added_fields.add(name)
                    if fasta_enabled:
                        attributes4 += f';{name}={value}'

    # Extract sequence if genomic file is provided
    gene_sequence = ""
    if seq_dict:
        if "gene" in include or "all" in include:
            gene_sequence = extract_sequence(seq_dict, chrom, *gene_positions.get(gene_id, [mRNA_start, mRNA_end]), strand)

    # Add gene entry only if it's a new gene
    if gene_id not in gene_ids_seen:
        gene_ids_seen.add(gene_id)
        start, end = gene_positions[gene_id]
        attributes_gene = f'gene_id "{gene_id}"; transcript_id ""; db_xref "{gene_id}"; gbkey "Gene"; gene_biotype "protein_coding"; '
        attributes_gene += attributes2
        attributes_gene += attributes3
        if str(gene_sequence) != "" and ("gene" in include or "all" in include):
            attributes_gene += f'sequence "{gene_sequence}"; '
        gene_line = f'{chrom}\tHANNO\tgene\t{start+1}\t{end}\t.\t{strand}\t.\t{attributes_gene}'
        gtf_lines.append(gene_line)
    
    # Initialize to collect exon and CDS sequences for this transcript
    concatenated_exon_sequence = ""
    complete_coding_sequence = ""

    # Add transcript entry
    attributes_transcript = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; experiment "-"; gbkey "mRNA"; model_evidence "HANNO"; transcript_biotype "mRNA"; '
    attributes_transcript += attributes2

    # Initialize ints for CDS frame calculation
    cumulative_cds_length = 0
    next_frame = 0

    # Reverse order if on - strand
    exon_indices = list(range(block_count))
    if strand == '-':
        exon_indices = list(reversed(exon_indices))

    exon_number = 1
    for i in exon_indices:
        exon_start = mRNA_start + block_starts[i] + 1
        exon_end = exon_start + block_len[i] - 1

        # Extract the sequence snippet for the exon if genomic file is provided
        exon_sequence = extract_sequence(seq_dict, chrom, exon_start - 1, exon_end, strand) if seq_dict else ""
        # Add this exon sequence to the concatenated sequence
        concatenated_exon_sequence += exon_sequence

        # Add exon entry
        attributes_exon = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; experiment "-";  model_evidence "HANNO"; transcript_biotype "mRNA"; exon_number "{exon_number}"; '
        attributes_exon += attributes2
        attributes_exon += attributes3
        if str(exon_sequence) != "" and ("exon" in include or "all" in include):
            attributes_exon += f'sequence "{exon_sequence}"; '
        exon_line = f'{chrom}\tHANNO\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{attributes_exon}'
        gtf_lines.append(exon_line)

        if exon_start <= cds_end and exon_end >= cds_start:
            cds_start_in_exon = max(exon_start, cds_start)
            cds_end_in_exon = min(exon_end, cds_end)

            # Extract sequence snippet for the CDS if genomic file is provided
            cds_sequence = extract_sequence(seq_dict, chrom, cds_start_in_exon - 1, cds_end_in_exon, strand) if seq_dict else ""
            complete_coding_sequence += cds_sequence

            # Add CDS entry
            frame = next_frame
            protein_id = transcript_id + '.p1'
            attributes_cds = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; gbkey "CDS"; protein_id "{protein_id}"; exon_number "{exon_number}"; '
            attributes_cds += attributes2
            attributes_cds += attributes3
            if str(cds_sequence) != "" and ("CDS" in include or "all" in include):
                attributes_cds += f'sequence "{cds_sequence}"; '
            cds_line = f'{chrom}\tHANNO\tCDS\t{cds_start_in_exon}\t{cds_end_in_exon}\t.\t{strand}\t{frame}\t{attributes_cds}'
            gtf_lines.append(cds_line)

            segment_length = cds_end_in_exon - cds_start_in_exon + 1
            cumulative_cds_length += segment_length

            next_frame = 3 - cumulative_cds_length % 3
            if next_frame == 3:
                next_frame = 0

        exon_number += 1

    if concatenated_exon_sequence:
        attributes_transcript += f'mRNA_sequence "{concatenated_exon_sequence}"; coding_sequence "{complete_coding_sequence}"; '
    # Now append the transcript line after processing the exons
    transcript_line = f'{chrom}\tHANNO\ttranscript\t{mRNA_start+1}\t{mRNA_end}\t.\t{strand}\t.\t{attributes_transcript}'
    gtf_lines.insert(1, transcript_line)  # Insert transcript line after gene entry

    # Write sequences to FASTA files if requested
    if fasta_enabled and input_file_path and mrna_fasta_path and cds_fasta_path:
        # Prepare the FASTA header
        header = f">{transcript_id} [gene={row.iloc[25]}] [product={row.iloc[14]}]"
        if attributes4 != "":
            header += f" [{attributes4.replace(';', ' ')}]"

        # Write mRNA sequence
        write_fasta_with_formatting(mrna_fasta_path, header, concatenated_exon_sequence)

        # Write CDS sequence
        write_fasta_with_formatting(cds_fasta_path, header, complete_coding_sequence)

    return "\n".join(gtf_lines)

# Function to convert one bedDB entry to GFF format
def bedDB_to_gff(row, field_map, column_names, gene_positions, gene_ids_seen, seq_dict=None, include=None, fasta_enabled=False, input_file_path=None, mrna_fasta_path=None, cds_fasta_path=None):
    gff_lines = []

    chrom = row.iloc[0]
    mRNA_start = int(row.iloc[1])
    mRNA_end = int(row.iloc[2])
    strand = row.iloc[5]
    cds_start = int(row.iloc[6] + 1)  # +1 due to GTF format being 1-indexed, unlike 0-indexed bed format
    cds_end = int(row.iloc[7])
    block_count = int(row.iloc[9])
    block_len = [int(x) for x in row.iloc[10].split(',') if x]
    block_starts = [int(x) for x in row.iloc[11].split(',') if x]

    # Check if exon sizes and starts match
    if len(block_len) != len(block_starts):
        raise ValueError("Problem with number of exons in bedDB")

    gene_id = row.iloc[3].rsplit('.', 1)[0]
    transcript_id = row.iloc[3]
    attributes2 = ""
    attributes3 = ""

    # Handle --add_fields and --add_field_names
    if not field_map:
        attributes2 += f';gene={row.iloc[25]};product={row.iloc[14]}'
    else:
        added_fields = set()
        for index, name in field_map.items():
            if index < len(row):
                value = row.iloc[index]
                if name not in added_fields:
                    attributes2 += f';{name}={value}'
                    added_fields.add(name)

    # Extract sequence if genomic file is provided
    gene_sequence = ""
    if seq_dict:
        if "gene" in include or "all" in include:
            gene_sequence = extract_sequence(seq_dict, chrom, *gene_positions.get(gene_id, [mRNA_start, mRNA_end]), strand)

    # Add gene entry only if it's a new gene
    if gene_id not in gene_ids_seen:
        gene_ids_seen.add(gene_id)
        start, end = gene_positions[gene_id]
        attributes_gene = f'ID=gene-{gene_id};Dbxref=GeneID:{gene_id};Name={gene_id};gbkey=Gene;gene_biotype=protein coding'
        attributes_gene += attributes2
        attributes_gene += attributes3
        if str(gene_sequence) != "" and ("gene" in include or "all" in include):
            attributes_gene += f';sequence={gene_sequence}'
        gene_line = f'{chrom}\tHANNO\tgene\t{start+1}\t{end}\t.\t{strand}\t.\t{attributes_gene}'
        gff_lines.append(gene_line)
    
    # Initialize to collect exon sequences for this transcript
    concatenated_exon_sequence = ""
    complete_coding_sequence = ""

    # Add transcript (mRNA) entry
    attributes_transcript = f'ID=mRNA-{transcript_id};Parent=gene-{gene_id};Dbxref=GeneID:{transcript_id};Name={transcript_id};gbkey=mRNA;model_evidence=HANNO;transcript_id={transcript_id}'
    attributes_transcript += attributes2

    # Initialize ints for CDS frame calculation
    cumulative_cds_length = 0
    next_frame = 0

    # Reverse order if on - strand
    exon_indices = list(range(block_count))
    if strand == '-':
        exon_indices = list(reversed(exon_indices))

    exon_number = 1
    for i in exon_indices:
        exon_start = mRNA_start + block_starts[i] + 1
        exon_end = exon_start + block_len[i] - 1

        # Extract the sequence snippet for the exon if genomic file is provided
        exon_sequence = extract_sequence(seq_dict, chrom, exon_start - 1, exon_end, strand) if seq_dict else ""
        concatenated_exon_sequence += exon_sequence # Add this exon sequence to the concatenated sequence

        # Add exon entry
        attributes_exon = f'ID=exon-{transcript_id}.{exon_number};Parent=mRNA-{transcript_id};Dbxref=GeneID:{gene_id};experiment=-;gbkey=mRNA;transcript_id={transcript_id}'
        attributes_exon += attributes2
        attributes_exon += attributes3
        if str(exon_sequence) != "" and ("exon" in include or "all" in include):
            attributes_exon += f';sequence={exon_sequence}'
        exon_line = f'{chrom}\tHANNO\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{attributes_exon}'
        gff_lines.append(exon_line)

        if exon_start <= cds_end and exon_end >= cds_start:
            cds_start_in_exon = max(exon_start, cds_start)
            cds_end_in_exon = min(exon_end, cds_end)

            # Extract sequence snippet for the CDS if genomic file is provided
            cds_sequence = extract_sequence(seq_dict, chrom, cds_start_in_exon - 1, cds_end_in_exon, strand) if seq_dict else ""
            complete_coding_sequence += cds_sequence # Add this CDS sequence to the concatenated sequence

            # Add CDS entry
            frame = next_frame
            protein_id = transcript_id + '.p1'
            attributes_cds = f'ID=cds-{transcript_id};Parent=mRNA-{transcript_id};Dbxref=GeneID:{gene_id};Name={gene_id};gbkey=CDS;protein_id={protein_id}'
            attributes_cds += attributes2
            attributes_cds += attributes3
            if str(cds_sequence) != "" and ("CDS" in include or "all" in include):
                attributes_cds += f';sequence={cds_sequence}'
            cds_line = f'{chrom}\tHANNO\tCDS\t{cds_start_in_exon}\t{cds_end_in_exon}\t.\t{strand}\t{frame}\t{attributes_cds}'
            gff_lines.append(cds_line)

            segment_length = cds_end_in_exon - cds_start_in_exon + 1
            cumulative_cds_length += segment_length

            next_frame = 3 - cumulative_cds_length % 3
            if next_frame == 3:
                next_frame = 0

        exon_number += 1

    # Add the concatenated exon sequence to the transcript attributes
    if concatenated_exon_sequence:
        attributes_transcript += f';mRNA_sequence={concatenated_exon_sequence};coding_sequence={complete_coding_sequence}'
    # Now append the transcript line after processing the exons
    transcript_line = f'{chrom}\tHANNO\tmRNA\t{mRNA_start+1}\t{mRNA_end}\t.\t{strand}\t.\t{attributes_transcript}'
    gff_lines.insert(1, transcript_line)  # Insert transcript line after gene entry

    # Write sequences to FASTA files if requested
    if fasta_enabled and input_file_path and mrna_fasta_path and cds_fasta_path:

        # Assuming you already have these variables defined in your function
        header = f">{transcript_id} [gene={row.iloc[25]}] [product={row.iloc[14]}]"
        if attributes3 != "":
            header+= f" [{attributes3.replace(';', ' ')}]"
    
        # Write mRNA sequence
        write_fasta_with_formatting(mrna_fasta_path, header, concatenated_exon_sequence)

        # Write CDS sequence
        write_fasta_with_formatting(cds_fasta_path, header, complete_coding_sequence)

    return "\n".join(gff_lines)

# Function to convert BED12 entry to GTF format
def bed12_to_gtf(row, field_map, seq_dict=None, include=None, fasta_enabled=False, input_file_path=None, mrna_fasta_path=None, cds_fasta_path=None):
    gtf_lines = []
    chrom = row.iloc[0]             # row.iloc[] was used to avoid problems in future pandas versions
    chrom_start = int(row.iloc[1])
    chrom_end = int(row.iloc[2])
    name = row.iloc[3]
    strand = row.iloc[5]
    thick_start = int(row.iloc[6] + 1)  # +1 due to GTF format being 1-indexed, unlike 0-indexed bed format
    thick_end = int(row.iloc[7])
    itemRGB = row.iloc[8]
    block_count = int(row.iloc[9])
    block_len = [int(x) for x in row.iloc[10].split(',') if x]
    block_starts = [int(x) for x in row.iloc[11].split(',') if x]

    # Check if exon sizes and starts match
    if len(block_len) != len(block_starts):
        raise ValueError("Problem with number of exons in bed12")

    # Create attributes
    gene_id = row.iloc[3].rsplit('.', 1)[0]
    transcript_id = name
    db_xref = name

    # Add itemRGB if provided or default
    added_fields = set()  # Track added fields to avoid duplicates
    attributes2 = ""
    attributes3 = ""
    if field_map:
        if 8 in field_map:
            attributes2 += f'{field_map[8]} "{itemRGB}"; '
            added_fields.add(field_map[8])
            if fasta_enabled:
                attributes3 += f';{field_map[8]}={itemRGB}'
        else:
            print("This is not a valid custom field entry for BED12, see --fields for available fields")
            sys.exit(0)

    # Add gene entry
    gene_sequence = ""
    if seq_dict:
        if "gene" in include or "all" in include:
            gene_sequence = extract_sequence(seq_dict, chrom, chrom_start, chrom_end, strand)
    attributes_gene = f'gene_id "{gene_id}"; transcript_id ""; db_xref "{gene_id}"; gbkey "Gene"; gene_biotype "protein_coding"; '
    attributes_gene += attributes2
    if str(gene_sequence) != "" and ("gene" in include or "all" in include):
            attributes_gene += f'sequence "{gene_sequence}"; '
    gene_line = f'{chrom}\tHANNO\tgene\t{chrom_start+1}\t{chrom_end}\t.\t{strand}\t.\t{attributes_gene}'
    gtf_lines.append(gene_line)

    # Initialize to collect exon and CDS sequences for this transcript
    concatenated_exon_sequence = ""
    complete_coding_sequence = ""

    # Add transcript entry
    attributes_transcript = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; experiment "-"; gbkey "mRNA";  model_evidence "HANNO"; transcript_biotype "mRNA"; '
    attributes_transcript += attributes2

    # Initialize ints for CDS frame calculation
    cumulative_cds_length = 0
    next_frame = 0

    exon_indices = list(range(block_count))
    if strand == '-':
        exon_indices = list(reversed(exon_indices))

    exon_number = 1
    for i in exon_indices:
        exon_start = chrom_start + block_starts[i] + 1
        exon_end = exon_start + block_len[i] - 1

        # Extract the sequence snippet for the exon if genomic file is provided
        exon_sequence = extract_sequence(seq_dict, chrom, exon_start - 1, exon_end, strand) if seq_dict else ""

        # Add this exon sequence to the concatenated sequence
        concatenated_exon_sequence += exon_sequence

        # Add exon entry
        attributes_exon = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; experiment "-";  model_evidence "HANNO"; transcript_biotype "mRNA"; exon_number "{exon_number}"; '
        attributes_exon += attributes2
        if str(exon_sequence) != "" and ("exon" in include or "all" in include):
            attributes_exon += f'sequence "{exon_sequence}"; '
        exon_line = f'{chrom}\tHANNO\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{attributes_exon}'
        gtf_lines.append(exon_line)

        if exon_start <= thick_end and exon_end >= thick_start:
            cds_start_in_exon = max(exon_start, thick_start)
            cds_end_in_exon = min(exon_end, thick_end)

            # Extract sequence snippet for the CDS if genomic file is provided
            cds_sequence = extract_sequence(seq_dict, chrom, cds_start_in_exon - 1, cds_end_in_exon, strand) if seq_dict else ""
            complete_coding_sequence += cds_sequence

            # Add CDS entry
            frame = next_frame
            protein_id = transcript_id + '.p1'
            attributes_CDS = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; db_xref "{db_xref}"; gbkey "CDS"; protein_id "{protein_id}"; exon_number "{exon_number}"; '
            attributes_CDS += attributes2
            if str(cds_sequence) != "" and ("CDS" in include or "all" in include):
                attributes_cds += f'sequence "{cds_sequence}"; '
            cds_line = f'{chrom}\tHANNO\tCDS\t{cds_start_in_exon}\t{cds_end_in_exon}\t.\t{strand}\t{frame}\t{attributes_CDS}'
            gtf_lines.append(cds_line)

            segment_length = cds_end_in_exon - cds_start_in_exon + 1
            cumulative_cds_length += segment_length

            next_frame = 3 - cumulative_cds_length % 3
            if next_frame == 3:
                next_frame = 0

        exon_number += 1

    if concatenated_exon_sequence:
        attributes_transcript += f'mRNA_sequence "{concatenated_exon_sequence}"; coding_sequence "{complete_coding_sequence}"; '
    # Now append the transcript line after processing the exons
    transcript_line = f'{chrom}\tHANNO\ttranscript\t{chrom_start+1}\t{chrom_end}\t.\t{strand}\t.\t{attributes_transcript}'
    gtf_lines.insert(1, transcript_line)  # Insert transcript line after gene entry

    # Write sequences to FASTA files if requested
    if fasta_enabled and input_file_path and mrna_fasta_path and cds_fasta_path:
        header = f">{transcript_id} [gene={gene_id}] [product=unknown]"
        if attributes3 != "":
            header += f" [{attributes3.replace(';', ' ')}]"
    
        # Write mRNA sequence
        write_fasta_with_formatting(mrna_fasta_path, header, concatenated_exon_sequence)

        # Write CDS sequence
        write_fasta_with_formatting(cds_fasta_path, header, complete_coding_sequence)

    return "\n".join(gtf_lines)

# Function to convert BED12 entry to GFF format
def bed12_to_gff(row, field_map, seq_dict=None, include=None, fasta_enabled=False, input_file_path=None, mrna_fasta_path=None, cds_fasta_path=None):
    gff_lines = []
    chrom = row.iloc[0]             # row.iloc[] was used to avoid problems in future pandas versions
    chrom_start = int(row.iloc[1])
    chrom_end = int(row.iloc[2])
    name = row.iloc[3]
    strand = row.iloc[5]
    thick_start = int(row.iloc[6] + 1)  # +1 due to GTF format being 1-indexed, unlike 0-indexed bed format
    thick_end = int(row.iloc[7])
    itemRGB = row.iloc[8]
    block_count = int(row.iloc[9])
    block_len = [int(x) for x in row.iloc[10].split(',') if x]
    block_starts = [int(x) for x in row.iloc[11].split(',') if x]

    # Check if exon sizes and starts match
    if len(block_len) != len(block_starts):
        raise ValueError("Problem with number of exons in bed12")

    # Create attributes
    gene_id = row.iloc[3].rsplit('.', 1)[0]
    transcript_id = name

    # Add itemRGB if provided or default to none
    added_fields = set()  # Track added fields to avoid duplicates
    attributes2=""
    if field_map:
        if 8 in field_map:
            attributes2 += f';{field_map[8]}={itemRGB}'
            added_fields.add(field_map[8])
        else:
            print("This is not a valid custom field entry for BED12, see --fields for available fields")
            sys.exit(0)

    # Extract sequence if genomic file is provided
    gene_sequence = ""
    if seq_dict:
        if "gene" in include or "all" in include:
            gene_sequence = extract_sequence(seq_dict, chrom, chrom_start, chrom_end, strand)

    # Add gene entry
    attributes_gene = f'ID=gene-{gene_id};Dbxref=GeneID:{gene_id};Name={gene_id};gbkey=Gene;gene_biotype=protein coding'
    attributes_gene += attributes2
    if str(gene_sequence) != "" and ("gene" in include or "all" in include):
            attributes_gene += f';sequence={gene_sequence}'
    gene_line = f'{chrom}\tHANNO\tgene\t{chrom_start+1}\t{chrom_end}\t.\t{strand}\t.\t{attributes_gene}'
    gff_lines.append(gene_line)

    # Initialize to collect exon sequences for this transcript
    concatenated_exon_sequence = ""
    complete_coding_sequence = ""

    # Add transcript entry
    attributes_transcript = f'ID=mRNA-{gene_id};Parent=gene-{gene_id};Dbxref=GeneID:{gene_id};Name={transcript_id};gbkey="mRNA";model_evidence=HANNO;transcript_id={transcript_id}'
    attributes_transcript += attributes2

    # Initialize ints for CDS frame calculation
    cumulative_cds_length = 0
    next_frame = 0

    exon_indices = list(range(block_count))
    if strand == '-':
        exon_indices = list(reversed(exon_indices))

    exon_number = 1
    for i in exon_indices:
        exon_start = chrom_start + block_starts[i] + 1
        exon_end = exon_start + block_len[i] - 1#

        # Extract the sequence snippet for the exon if genomic file is provided
        exon_sequence = extract_sequence(seq_dict, chrom, exon_start - 1, exon_end, strand) if seq_dict else ""
        concatenated_exon_sequence += exon_sequence # Add this exon sequence to the concatenated sequence

        # Add exon entry
        attributes_exon = f'ID=exon-{transcript_id}.{exon_number};Parent=mRNA-{transcript_id};Dbxref=GeneID:{gene_id};experiment=-;gbkey=mRNA'
        attributes_exon += attributes2
        if str(exon_sequence) != "" and ("exon" in include or "all" in include):
            attributes_exon += f';sequence={exon_sequence}'
        exon_line = f'{chrom}\tHANNO\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{attributes_exon}'
        gff_lines.append(exon_line)

        if exon_start <= thick_end and exon_end >= thick_start:
            cds_start_in_exon = max(exon_start, thick_start)
            cds_end_in_exon = min(exon_end, thick_end)
            
            # Extract sequence snippet for the CDS if genomic file is provided
            cds_sequence = extract_sequence(seq_dict, chrom, cds_start_in_exon - 1, cds_end_in_exon, strand) if seq_dict else ""
            complete_coding_sequence += cds_sequence # Add this CDS sequence to the concatenated sequence

            # Add CDS entry
            frame = next_frame
            attributes_CDS = f'ID=cds-{transcript_id};Parent=mRNA-{transcript_id};Dbxref=GeneID:{gene_id};Name={gene_id};gbkey=CDS'
            attributes_CDS += attributes2
            if str(cds_sequence) != "" and ("CDS" in include or "all" in include):
                attributes_cds += f';sequence={cds_sequence}'
            cds_line = f'{chrom}\tHANNO\tCDS\t{cds_start_in_exon}\t{cds_end_in_exon}\t.\t{strand}\t{frame}\t{attributes_CDS}'
            gff_lines.append(cds_line)

            segment_length = cds_end_in_exon - cds_start_in_exon + 1
            cumulative_cds_length += segment_length

            next_frame = 3 - cumulative_cds_length % 3
            if next_frame == 3:
                next_frame = 0

        exon_number += 1
        
    # Add the concatenated exon sequence to the transcript attributes
    if concatenated_exon_sequence:
        attributes_transcript += f';mRNA_sequence={concatenated_exon_sequence};coding_sequence={complete_coding_sequence}'
    # Now append the transcript line after processing the exons
    transcript_line = f'{chrom}\tHANNO\tmRNA\t{chrom_start+1}\t{chrom_end}\t.\t{strand}\t.\t{attributes_transcript}'
    gff_lines.insert(1, transcript_line)  # Insert transcript line after gene entry

    # Write sequences to FASTA files if requested
    if fasta_enabled and input_file_path and mrna_fasta_path and cds_fasta_path:
        header = f">{transcript_id} [gene={gene_id}] [product=unknown]"
        if attributes2 != "":
            header += f" [{attributes2.replace(';', ' ')}]"
    
        # Write mRNA sequence
        write_fasta_with_formatting(mrna_fasta_path, header, concatenated_exon_sequence)

        # Write CDS sequence
        write_fasta_with_formatting(cds_fasta_path, header, complete_coding_sequence)

    return "\n".join(gff_lines)

# Function to convert BED6 entry to GTF/GFF format
def bed6_to_gtf(row, mode):
    chrom = row.iloc[0]
    start = int(row.iloc[1]) + 1  # BED format is 0-based, GTF is 1-based
    end = int(row.iloc[2])
    name = row.iloc[3]
    score = row.iloc[4]
    strand = row.iloc[5]
    
    if mode == 'gff':
        attributes = f'ID={name};Name={name}'
    else:
        attributes = f'gene_id "{name}"; transcript_id "{name}";'
    
    transcript_line = f'{chrom}\tBED\tregion\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}'
    return transcript_line

# Function to convert BED5 entry to GTF/GFF format
def bed5_to_gtf(row, mode):
    chrom = row.iloc[0]
    start = int(row.iloc[1]) + 1
    end = int(row.iloc[2])
    name = row.iloc[3]
    score = row.iloc[4]
    
    if mode == 'gff':
        attributes = f'ID={name};Name={name}'
    else:
        attributes = f'gene_id "{name}"; transcript_id "{name}";'
    
    gtf_line = f'{chrom}\tBED\tregion\t{start}\t{end}\t{score}\t.\t.\t{attributes}'
    return gtf_line

# Function to convert BED4 entry to GTF/GFF format
def bed4_to_gtf(row, mode):
    chrom = row.iloc[0]
    start = int(row.iloc[1]) + 1
    end = int(row.iloc[2])
    name = row.iloc[3]
    
    if mode == 'gff':
        attributes = f'ID={name};Name={name}'
    else:
        attributes = f'gene_id "{name}"; transcript_id "{name}";'
    
    gtf_line = f'{chrom}\tBED\tregion\t{start}\t{end}\t.\t.\t.\t{attributes}'
    return gtf_line

# Helper function to check file type based on column count and data types
def guess_file_type(df):
    col_count = df.shape[1]
    
    if col_count == 4:
        if df.dtypes.iloc[0] =='object' and df.dtypes.iloc[1] == 'int64' and df.dtypes.iloc[2] == 'int64' and df.dtypes.iloc[3] == 'object':
            return 'bed4'
    elif col_count == 5:
        if df.dtypes.iloc[0] =='object' and df.dtypes.iloc[1] == 'int64' and df.dtypes.iloc[2] == 'int64' and df.dtypes.iloc[3] == 'object':
            return 'bed5'
    elif col_count == 6:
        if df.dtypes.iloc[0] =='object' and df.dtypes.iloc[1] == 'int64' and df.dtypes.iloc[2] == 'int64'and df.dtypes.iloc[3] == 'object':
            strand_values = df.iloc[:, 5].unique()
            if set(strand_values).issubset({'+', '-'}):
                return 'bed6'
    elif col_count == 12:
        # Check if it's likely BED12 format by validating key columns
        if df.dtypes.iloc[1] == 'int64' and df.dtypes.iloc[2] == 'int64' and df.dtypes.iloc[5] == 'object':
            strand_values = df.iloc[:, 5].unique()
            if set(strand_values).issubset({'+', '-'}):
                return 'bed12'
    elif col_count == 46:
        # Check if it's likely bedDB format
        if df.dtypes.iloc[6] == 'int64' and df.dtypes.iloc[5] == 'object':
            strand_values = df.iloc[:, 5].unique()
            if set(strand_values).issubset({'+', '-'}):
                return 'bedDB'
    
    return 'unknown'

# Main function
def main():
    # Setup argparse for command-line options
    parser = argparse.ArgumentParser(description='Convert bedDB or BED4/5/6/12 file file to GTF/GFF format.', add_help=False)
    parser.add_argument('input_file', type=str, nargs='?', help='Input bedDB or BED4/5/6/12 file with complete path')
    parser.add_argument('-o', '--output', type=str, help='Output GTF file name (optional, auto-generated from input file name if not provided)')
    parser.add_argument('-h', '--help', action='store_true', help='Show help message')
    parser.add_argument('--add_fields', type=str, help='Comma-separated list of field indices to include in GTF attributes, default is 26,15')
    parser.add_argument('--add_field_names', type=str, help='Comma-separated list of custom names for the fields in GTF attributes')
    parser.add_argument('--fields', action='store_true', help='Shows all additional fields from bedDB or BED12 that can be added to GTF attributes with index, an input file must be provided')
    parser.add_argument('--gff', action='store_true', help='The script produces a GFF file, instead of a GTF file' )
    parser.add_argument('--genome', type=str, help='Optional second input file with complete path (genomic sequences)')
    parser.add_argument('--fasta', action='store_true', help='Optional: Generate two FASTA files, one for mRNA and one for CDS sequences, with default filenames.')
    # Add dependent argument
    parser.add_argument('--include', type=str, nargs='+', choices=['gene', 'CDS', 'exon', 'all'],
                        help='Specify which entries to include sequences for: gene, CDS, exon, or all. If not provided, sequences will not be added. Use space to separate multiple options.')
    
    # Parse arguments
    args = parser.parse_args()

    # Display help if '-h' or '--help' is passed
    if args.help:
        print("Program function: Convert bedDB or BED4/5/6/12 file to GTF/GFF format")
        print("Arguments:")
        print("  input_file         Input bedDB or BED4/5/6/12 file with complete path")
        print("  -o, --output       Output GTF file name (optional, auto-generated from input file name if not provided)")
        print("  --add_fields       Comma-separated list of field indices to include in attributes, default is 26,15.")
        print("  --add_field_names  Comma-separated list of custom names for the fields in attributes")
        print("""  --fields           Show all additional fields from bedDB or BED12 that can be added to GTF attributes with index
                     (input bedDB or BED12 file necessary)""")
        print("  --gff              The script will produce a GFF file, instead of a GTF file")
        print("  --genome           Optional input genomic file(fasta format) with complete path")
        print("""  --include          Specify which entries to include sequences for: gene, CDS, exon, or all. If not provided, 
                      sequences will not be added. Use space to separate multiple options.""")
        print("  --fasta             Optional: Generate two FASTA files, one for mRNA and one for CDS sequences, with default filenames.")
        print("Incase of 'ModuleNotFoundError: No module named 'pandas'', use 'pip install pandas' to get the module")
        print("  -h, --help         Show this help message and exit")
        sys.exit(0)

    column_names=[]

    # Ensure an input file is provided
    if not args.input_file:
        print("Error: No input file provided.")
        print("Use -h or --help for usage information.")
        sys.exit(1)

    # Automatically generate output file name if not provided
    input_dir = os.path.dirname(os.path.abspath(args.input_file))  # Get the directory of the input file
    output_file= ""
    if not args.output:
        input_basename = os.path.splitext(args.input_file)[0]
        field_type = 'custom' if args.add_fields else 'long'  # Determine field type based on user selection
        seq_status = 'with_seq' if args.genome else 'no_seq'  # Determine sequence inclusion

        # Generate the file name based on field type, seq status, and file format (GTF or GFF)
        if args.gff:
            output_file = f"{input_basename}-{field_type}-{seq_status}.gff"
        else:
            output_file = f"{input_basename}-{field_type}-{seq_status}.gtf"
    else:
        # Use the provided output file name, but ensure it is written to the input directory
        output_basename = os.path.basename(args.output)  # Get the file name part of the output
        output_file = os.path.join(input_dir, output_basename)  # Combine input directory with the provided output file name

    # Initialize FASTA files if --fasta option is enabled
    mrna_fasta_path = None
    cds_fasta_path = None
    if args.fasta:
        mrna_fasta_path, cds_fasta_path = initialize_fasta_files(input_dir, output_file)

    # Handle custom field additions
    field_map = {}
    if args.add_fields and args.add_field_names:
        field_indices = [int(idx) - 1 for idx in args.add_fields.split(',')]  # 0-based indexing
        field_names = args.add_field_names.split(',')
        if len(field_indices) != len(field_names):
            sys.exit("Error: The number of field indices and field names must be the same.")
        field_map = dict(zip(field_indices, field_names))
    elif args.add_fields or args.add_field_names:
        sys.exit("Error: --add_fields and --add_field_names must be used together.")

    seq_dict = {}
    if args.genome:
        seq_dict = load_genomic_sequences(args.genome)

    mode = 'gtf'
    if args.gff:
        mode='gff'

    # Try reading the file first without a header, and then with a header if necessary
    try:
        # First, attempt to read without a header
        bed_df = pd.read_csv(args.input_file, sep='\t', header=None, low_memory=False)
        file_type = guess_file_type(bed_df)
        
        if file_type == 'unknown':
            # If not recognized as valid, try reading with a header (likely for bedDB)
            bed_df = pd.read_csv(args.input_file, sep='\t', header='infer')
            file_type = guess_file_type(bed_df)
            column_names = bed_df.columns.tolist()

        if file_type == 'bed4':
            print(f"Detected bed4 file: {args.input_file}")
            lines = bed_df.apply(lambda row: bed4_to_gtf(row, mode), axis=1)
        elif file_type == 'bed5':
            print(f"Detected bed5 file: {args.input_file}")
            lines = bed_df.apply(lambda row: bed5_to_gtf(row, mode), axis=1)
        elif file_type == 'bed6':
            print(f"Detected bed6 file: {args.input_file}")
            lines = bed_df.apply(lambda row: bed6_to_gtf(row, mode), axis=1)
        elif file_type == 'bed12':
            print(f"Detected bed12 file: {args.input_file}")
            if args.fields:
                print("9: itemRGB")
                sys.exit(0)
            else:
                converter = bed12_to_gff if args.gff else bed12_to_gtf
                include = args.include if args.include else []
                lines = bed_df.apply(lambda row: converter(row, field_map, seq_dict, include, args.fasta, args.input_file, mrna_fasta_path, cds_fasta_path), axis=1)
        elif file_type == 'bedDB':
            print(f"Detected bedDB file: {args.input_file}")
            if args.fields:
                if not column_names:
                    # If column names haven't been populated yet, read the file with header again
                    bed_df = pd.read_csv(args.input_file, sep='\t', header='infer')
                    column_names = bed_df.columns.tolist()

                # Display all available fields with their index (1-based)
                print("Available fields in bedDB:")
                for i, col in enumerate(column_names[12:], start=13):
                    print(f"{i}: {col}")
                sys.exit(0)
            else:
                # Read file twice, to get correct gene coordinates
                # First pass: collect gene positions
                gene_positions = collect_gene_positions(bed_df)

                # Second pass: generate GTF lines
                gene_ids_seen = set()  # Keep track of processed gene IDs
                converter = bedDB_to_gff if args.gff else bedDB_to_gtf
                include = args.include if args.include else []
                lines = bed_df.apply(lambda row: converter(row, field_map, column_names, gene_positions, gene_ids_seen, seq_dict, include, args.fasta, args.input_file, mrna_fasta_path, cds_fasta_path), axis=1)
        else:
            print(f"Error: {args.input_file} does not appear to be a valid BED4, BED5, BED6, BED12 or bedDB file.")
            sys.exit(1)

    except FileNotFoundError:
        print(f"Error: File {args.input_file} not found.")
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Error: Could not parse {args.input_file} as a valid tab-separated file.")
        sys.exit(1)

    # Write output file
    with open(output_file, "w", newline='') as f:
        if args.gff:
            f.write(f"##gff-version 3\n")
        for line in lines:
            f.write(f"{line}\n")


    print(f"File created: {output_file}")

if __name__ == '__main__':
    main()
