##########
### 07 january 2025
### leanie 
### operonandon start
###
### daRulez:
### 1. same strand
### 2. serially located
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
    
def find_operons(df, sep_thresh):
    """ 
    Identify operons based on daRulez:
    1. co-orientation 
    2. seriality of features. 

    --sep_thresh 
    default = 500 (0-500)
    separation between features
    how the user defines inter-genic regions
    """
    operons = []
    df['attributes_dict'] = df['attribute'].apply(parse_attributes)
    
    sorted_df = df.sort_values(by=['seqname', 'start'])
    current_operon = []
    last_end = None
    last_strand = None
    
    for index, row in sorted_df.iterrows():
        if (current_operon and (row['strand'] != last_strand or (row['start'] - last_end > sep_thresh))):
            operons.append(current_operon)
            current_operon = []
        current_operon.append(row)
        last_end = row['end']
        last_strand = row['strand']

    if current_operon:
        operons.append(current_operon)

    return operons
    

def output_intergenic_regions(operons, sep_thresh, intergenic_bed_file):
    """ 
    Write intergenic regions to a BED file. 
    """
    with open(intergenic_bed_file, 'w') as f:
        last_end = None
        last_chrom = None
        last_operon_id = None ### keep track of the operon id for naming the region

        for operon in operons:
            current_operon_id = f"{operon[0]['attributes_dict'].get('Name', operon[0]['attributes_dict'].get('locus_tag', 'unknown'))}_Op"
            chrom = operon[0]['seqname']
            strand = operon[0]['strand']
            for feature in operon:
                chrom = feature['seqname']
                start = min(feature['start'] for feature in operon)#feature['start']
                if last_end is not None and last_chrom == chrom and start > last_end + 1:
                    intergenic_start = last_end + 1
                    intergenic_end = start - 1

                    region_name = f"{last_operon_id}_{current_operon_id}"
                    f.write(f"{chrom}\t{intergenic_start}\t{intergenic_end}\t{region_name}\t{strand}\n")
                
                last_end = max(feature['end'] for feature in operon)
                last_chrom = chrom
                last_operon_id = current_operon_id

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

def operons_to_csv(operons, sep_thresh, output_file):
    """ 
    Create operon results and return in CSV format.
    oriented: same strand (0 for -, 1 for +)
    oriented_nearby: same strand AND w/i sep_thresh (TRUE, 1)
    """
    written = set()
    with open(output_file, 'w') as f:
        ### header
        f.write("Feature, Locus_tag, Gene_type, Operon_ID, Operon_Start, Operon_End, Feature_Start, Feature_End, Strand, Oriented, Oriented_Nearby\n" )
           
        for operon in operons:
            
            operon_id = f"{operon[0]['attributes_dict'].get('Name', operon[0]['attributes_dict'].get('locus_tag', 'unknown'))}_Op"
            operon_start = min(feature['start'] for feature in operon)
            operon_end = max(feature['end'] for feature in operon)
            strand = operon[0]['strand']  # All features in operon have the same strand
            
            previous_end = None
            previous_strand = None
            
            f.write(f"Operon,{operon_id},{operon_id},,Operon,{operon_start},{operon_end},,,{strand},,\n")

            for feature in operon:
                
                attributes = feature['attributes_dict']
                locus_tag = attributes.get('locus_tag')
                gene_type = 'pseudogene' if 'pseudo' in attributes else 'gene'
                feature_name = attributes.get('Name')
                unique_key = (locus_tag, operon_id, strand)

                feature_start = feature['start']
                feature_end = feature['end']

                ### oriented and oriented_nearby
                oriented = 1 if strand == '+' else 0
                if previous_strand == strand and previous_end is not None and (feature['start'] - previous_end <= sep_thresh):
                    oriented_nearby = 'TRUE'
                else:
                    oriented_nearby = 'FALSE'

                ### write results to csv
                ### no duplicates
                if unique_key not in written:
                    f.write(f"{feature_name},{locus_tag}, {gene_type} ,{operon_id},{operon_start}, {operon_end}, {feature_start},{feature_end},{strand}, {oriented}, {oriented_nearby}\n")
                    written.add(unique_key)

                ### Update previous_end and previous_strand for the next iteration
                previous_end = feature['end']
                previous_strand = strand


def main():
    parser = argparse.ArgumentParser(description='Identify operons in a GFF file and output to CSV, filtered to genes and pseudogenes.')
    parser.add_argument('--gff_file', type=str, help='Path to the GFF file')
    parser.add_argument('--output_file', type=str, help='Output base_name for files')
    parser.add_argument('-b','--bed', action='store_true', help='Output operons in BED format' )
    parser.add_argument('--sep_thresh', type=int, default=500, choices=range(0, 501), help='Maximum number of bases between features to consider them as part of the same operon (0-500, default 500)')
    parser.add_argument('--inter', action='store_true', help='Output BED file for intergenic regions')
    #parser.add_argument('v', '--verbose', action='store_true')
    #parser.add_argument('-h', '--help', action='store_true')

    ### parse arguments
    args = parser.parse_args()
    
    ### name output and files
    base_name = os.path.splitext(args.output_file)[0]
    csv_file = f"{base_name}.csv"
    bed_file = f"{base_name}.bed"
    intergenic_bed_file = f"{base_name}_intergenic.bed"

    ### read gff 
    df = read_gff(args.gff_file)

    ### finding out operons
    ### taking user input 
    operons = find_operons(df, args.sep_thresh)
    
    ### output csv
    operons_to_csv(operons, args.sep_thresh, csv_file)

    ### for bed flag
    if args.bed:
        operons_2_bed(operons, bed_file)


    if args.inter:
        output_intergenic_regions(operons, args.sep_thresh, intergenic_bed_file)


if __name__ == '__main__':
    main()
