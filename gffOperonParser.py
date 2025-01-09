##########
### 07 january 2025
### leanie 
### operonandon start
###
### same strand
### serially located --> WIP: overlap <500 bp
### input: gff
### output: csv ??
##########
### use
### python3 script/operon.start.py <input.gff> <output.csv>
##########
##########
##########
##########

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

def find_operons(df):
    """ Identify operons based on daRulez: co-orientation and seriality of features. """
    operons = []
    df['attributes_dict'] = df['attribute'].apply(parse_attributes)
    
    sorted_df = df.sort_values(by=['seqname', 'start'])
    current_operon = []
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
                unique_key = (locus_tag, operon_id, strand)
                if unique_key not in written:
                    f.write(f"{feature_name},{locus_tag}, {gene_type} ,{operon_id},{operon_start}, {operon_end}, {strand}\n")
                    written.add(unique_key)


def main():
    parser = argparse.ArgumentParser(description='Identify operons in a GFF file and output to CSV, filtered to genes and pseudogenes.')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')
    
    args = parser.parse_args()
    
    df = read_gff(args.gff_file)
    operons = find_operons(df)
    operons_to_csv(operons, args.output_file)

if __name__ == '__main__':
    main()
