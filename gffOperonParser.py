##########
### 07 january 2025
### leanie 
### operonandon start
###
### daRulez:
### 1. same strand
### 2. serially located
###	- overlap of less than 500 --> WIP
### 3. no fight club
### 
### 
### input: gff
### output: csv +
##########
### use
### python3 script/gffOperonParser.py <input.gff> <output.csv>
##########
##########
##########
##########
import os
import pandas as pd
import argparse

def read_gff(gff_file):
    """ Read GFF file and return a DataFrame, ignoring genome """
    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    ### Filter out entire genome if feature type is region and covers the full genome length
    df = df[~((df['feature'] == 'region') & (df['start'] == 1) & (df['end'] == df['end'].max()))]
    ### filter to genes and pseudogenes. 
    df = df[df['feature'].isin(['gene', 'pseudogene'])]

    return df

def parse_attributes(attribute_str):
    """ Parse the attribute column and convert to a dict. """
    attributes = {}
    for item in attribute_str.split(';'):
        key, value = item.split('=')
        attributes[key] = value
    return attributes

def extract_gene_name(attributes):
    """ Extract gene name or gene ID from the gff attributes"""
    for attr in attributes.split(';'):
        if attr.startswith('Name'):
            return attr.split('=')[1]
        if attr.startswith('ID='):
            return attr.split('=')[1]
    return 'unknown' ## default for no name/id
    
def find_operons(df):
    """ Identify operons based on daRulez: co-orientation and seriality of features. """
    operons = []
    df['attributes_dict'] = df['attribute'].apply(parse_attributes)
    
    sorted_df = df.sort_values(by=['seqname', 'start'])
    current_operon = []
    last_end = None
    last_strand = None
    
    for index, row in sorted_df.iterrows():
        if (current_operon and (row['strand'] != last_strand)):
                operons.append(current_operon)
                current_operon = []
        current_operon.append(row)
        last_strand = row['strand']

        if current_operon:
            operons.append(current_operon)

    return operons
    
def gff_2_bed(gff_file, bed_file):
    with open(gff_file, 'r') as gff, open(bed_file, 'w') as bed:
        for line in gff:
            if line.startswith("#"):
                continue #skip header
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue #skip incomplete lines

            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            ### adjust for bed format starting at 0
            start_bed = int(start) - 1

            gene_name = extract_gene_name(attributes)

            bed_line = f"{chrom}\t{start_bed}\t{end}\t{gene_name}\t{score}\t{strand}\n"
            bed.write(bed_line)

def operons_2_bed(operons, bed_file):
    """ write operons to a bed file"""
    with open(bed_file, 'w') as bed:
        for operon in operons:
            chrom = operon[0]['seqname']
            start = min(feature['start'] for feature in operon) -1 ### 0 start
            end = max(feature['end'] for feature in operon)
            strand = operon[0]['strand']
            operon_id = f"{operon[0]['attributes_dict'].get('Name', operon[0]['attributes_dict'].get('locus_tag', 'unknown'))}_Op"
            score = 0 ### place holder 

            bed_line = f"{chrom}\t{start}\t{end}\t{operon_id}\t{score}\t{strand}\n"
            bed.write(bed_line)

def operons_to_csv(operons, output_file):
    """ Convert list of operons to CSV format. """
    written = set()
    with open(output_file, 'w') as f:
        ### header
        f.write("Feature, Locus_tag, Gene_type, Operon_ID, Operon_Start, Operon_End, Strand\n" )
           
        for operon in operons:
            operon_id = f"{operon[0]['attributes_dict'].get('Name', operon[0]['attributes_dict'].get('locus_tag', 'unknown'))}_Op"
            operon_start = min(feature['start'] for feature in operon)
            operon_end = max(feature['end'] for feature in operon)
            strand = operon[0]['strand']  # All features in operon have the same strand
            
            for feature in operon:
                attributes = feature['attributes_dict']
                locus_tag = attributes.get('locus_tag')
                gene_type = 'pseudogene' if 'pseudo' in attributes else 'gene'
                feature_name = attributes.get('Name')
                f.write(f"{feature_name},{locus_tag}, {gene_type} ,{operon_id},{operon_start}, {operon_end}, {strand}\n")

def main():
    parser = argparse.ArgumentParser(description='Identify operons in a GFF file and output to CSV, filtered to genes and pseudogenes.')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file')
    parser.add_argument('output_file', type=str, help='Output base_name for files')
    parser.add_argument('-b','--bed', action='store_true', help='Output operons in BED format' )
    #parser.add_argument('v', '--verbose', action='store_true')
    #parser.add_argument('-h', '--help', action='store_true')


    args = parser.parse_args()
    
    base_name = os.path.splitext(args.output_file)[0]
    csv_file = f"{base_name}.csv"
    bed_file = f"{base_name}.bed"

    df = read_gff(args.gff_file)
    operons = find_operons(df)
    
    ### output csv
    operons_to_csv(operons, csv_file)

    ### for bed flag
    if args.bed:
        operons_2_bed(operons, bed_file)


if __name__ == '__main__':
    main()
