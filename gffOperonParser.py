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
    
    ### seq+ strand
    for (seqname, strand), group in df.groupby(['seqname', 'strand']):
        sorted_group = group.sort_values(by = 'start')
        current_operon = []
        last_end = None

        for index, row in sorted_group.iterrows():
            if current_operon and (row['start'] - last_end > 500):
                operons.append(current_operon)
                current_operon = []
            current_operon.append(row)
            last_end = row['end']

        if current_operon:
            operons.append(current_operon)

    return operons
    


def operons_to_csv(operons, output_file):
    """ Convert list of operons to CSV format. """
    with open(output_file, 'w') as f:
        for operon in operons:
            attributes = operon[0]['attributes_dict']
            identifier = attributes.get('Name', attributes.get('locus_tag', 'unknown'))
            operon_id = identifier + '_Op'
            for feature in operon:
                attributes = feature['attributes_dict']
                f.write(f"{feature['seqname']},{feature['start']},{feature['end']},"
                        f"{attributes.get('gene', attributes.get('Name',''))},{feature['score']},{feature['strand']},{feature['feature']},"
                        f"{operon_id},"
                        f"{attributes.get('locus_tag','')},{attributes.get('product','')}\n")

def main():
    parser = argparse.ArgumentParser(description='Identify operons in a GFF file. Output results to CSV.')
    parser.add_argument('gff_file', type=str, help='Path to the GFF file')
    parser.add_argument('output_file', type=str, help='Path and name of results. Output CSV file')
    
    args = parser.parse_args()
    
    df = read_gff(args.gff_file)
    operons = find_operons(df)
    operons_to_csv(operons, args.output_file)

if __name__ == '__main__':
    main()
