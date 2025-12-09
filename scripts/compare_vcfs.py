import argparse
import os
import numpy as np
import pandas as pd

def indel_chunks_match(a, b):
    
    if a is None or b is None:
        return False
    
    a = str(a).upper()
    b = str(b).upper()
    
    if not a or not b:
        return False
    
    if a == b:
        return True
    if a in b or b in a:
        return True
    
    return False

def insertions_match(vcf1_ref, vcf1_alt, vcf2_ref, vcf2_alt):

    #prefixmatch
    v1_pre = vcf1_alt[len(vcf1_ref):]
    v2_pre = vcf2_alt[len(vcf2_ref):]

    #suffixmatch
    v1_suf = vcf1_alt[:-len(vcf1_ref)]
    v2_suf = vcf2_alt[:-len(vcf2_ref)]

    #match prefix with suffix
    return (
        indel_chunks_match(v1_pre, v2_pre) or
        indel_chunks_match(v1_suf, v2_suf) or
        indel_chunks_match(v1_pre, v2_suf) or
        indel_chunks_match(v1_suf, v2_pre)
    )

def deletions_match(vcf1_ref, vcf1_alt, vcf2_ref, vcf2_alt):
    
    #prefixmatch
    v1_pre = vcf1_ref[len(vcf1_alt):]
    v2_pre = vcf2_ref[len(vcf2_alt):]

    #suffixmatch
    v1_suf = vcf1_ref[:-len(vcf1_alt)]
    v2_suf = vcf2_ref[:-len(vcf2_alt)]

    
    #match prefix with suffix
    return (
        indel_chunks_match(v1_pre, v2_pre) or
        indel_chunks_match(v1_suf, v2_suf) or
        indel_chunks_match(v1_pre, v2_suf) or
        indel_chunks_match(v1_suf, v2_pre)
    )

def with_truth(truth, bcf_vcf, snippy_vcf):
    output = pd.DataFrame(columns=['pos','truth_ref','truth_alt','bcf_ref','bcf_alt','snippy_ref','snippy_alt'])
    window = 10 #max insertion/deletion length
    bcf_found = []
    snippy_found = []

    for i in range(len(truth)):
        pos_value = truth.at[i, 'Pos']
        truth_ref = truth.at[i, 'original']
        if pd.isna(truth_ref):
            truth_ref = "" 
        truth_alt = truth.at[i, 'new']
        if pd.isna(truth_alt):
            truth_alt = "" 

        bcf_row, bcf_ref, bcf_alt = None, None, None
        snippy_row, snippy_ref, snippy_alt = None, None, None

        #if snp
        if len(truth_ref) == 1 and len(truth_alt) == 1:
            #check bcf_tools 
            bcf_row = bcf_vcf.loc[bcf_vcf['Pos'] == pos_value]
            if not bcf_row.empty:
                bcf_ref = bcf_row['Ref'].iloc[0]
                bcf_alt = bcf_row['Alt'].iloc[0]
            else:
                bcf_ref = None
                bcf_alt = None
            #check snippy
            snippy_row = snippy_vcf.loc[snippy_vcf['Pos'] == pos_value]
            if not snippy_row.empty:
                snippy_ref = snippy_row['Ref'].iloc[0]
                snippy_alt = snippy_row['Alt'].iloc[0]
            else:
                snippy_ref = None
                snippy_alt = None

        #insertion
        if len(truth_ref) < len(truth_alt):
            min_pos = pos_value - window
            max_pos = pos_value + window
            #check bcf
            bcf_row = bcf_vcf.loc[(bcf_vcf['Pos'] >= min_pos) & (bcf_vcf['Pos'] <= max_pos)]
            if not bcf_row.empty:
                bcf_ref = bcf_row['Ref'].iloc[0]
                bcf_alt = bcf_row['Alt'].iloc[0]
                if len(str(bcf_alt)) <= len(str(bcf_ref)):
                    bcf_ref = None
                    bcf_alt = None
                else:
                    if not insertions_match(truth_ref, truth_alt, bcf_ref, bcf_alt):
                        bcf_ref = None
                        bcf_alt = None
            else:
                bcf_ref = None
                bcf_alt = None

            #check snippy
            snippy_row = snippy_vcf.loc[(snippy_vcf['Pos'] >= min_pos) & (snippy_vcf['Pos'] <= max_pos)]
            if not snippy_row.empty:
                snippy_ref = snippy_row['Ref'].iloc[0]
                snippy_alt = snippy_row['Alt'].iloc[0]
                if len(str(snippy_alt)) <= len(str(snippy_ref)):
                    snippy_ref = None
                    snippy_alt = None
                else:
                    if not insertions_match(truth_ref, truth_alt, snippy_ref, snippy_alt):
                        snippy_ref = None
                        snippy_alt = None
            else:
                snippy_ref = None
                snippy_alt = None


        #deletion
        if len(truth_ref) > len(truth_alt):
            min_pos = pos_value - window
            max_pos = pos_value + window
            #check bcf
            bcf_row = bcf_vcf.loc[(bcf_vcf['Pos'] >= min_pos) & (bcf_vcf['Pos'] <= max_pos)]
            if not bcf_row.empty:
                bcf_ref = bcf_row['Ref'].iloc[0]
                bcf_alt = bcf_row['Alt'].iloc[0]
                if len(str(bcf_alt)) >= len(str(bcf_ref)):
                    bcf_ref = None
                    bcf_alt = None
                else:
                    if not deletions_match(truth_ref, truth_alt, bcf_ref, bcf_alt):
                        bcf_ref = None
                        bcf_alt = None
            else:
                bcf_ref = None
                bcf_alt = None

            #check snippy
            snippy_row = snippy_vcf.loc[(snippy_vcf['Pos'] >= min_pos) & (snippy_vcf['Pos'] <= max_pos)]
            if not snippy_row.empty:
                snippy_ref = snippy_row['Ref'].iloc[0]
                snippy_alt = snippy_row['Alt'].iloc[0]
                if len(str(snippy_alt)) >= len(str(snippy_ref)):
                    snippy_ref = None
                    snippy_alt = None
                else:
                    if not deletions_match(truth_ref, truth_alt, snippy_ref, snippy_alt):
                        snippy_ref = None
                        snippy_alt = None
            else:
                snippy_ref = None
                snippy_alt = None

        if not (bcf_ref == None and bcf_alt == None):
            bcf_found.append(bcf_row['Pos'].iloc[0])

        if not(snippy_ref == None and snippy_alt == None):
            snippy_found.append(snippy_row['Pos'].iloc[0])

        output.loc[len(output)] = {
        'pos': pos_value,
        'truth_ref': truth_ref,
        'truth_alt': truth_alt,
        'bcf_ref': bcf_ref,
        'bcf_alt': bcf_alt,
        'snippy_ref': snippy_ref,
        'snippy_alt': snippy_alt}

    #bcf-unique calls
    for i in range(len(bcf_vcf)):
        row = bcf_vcf.iloc[i]
        row_pos = row['Pos']
        if row_pos not in bcf_found:
            output.loc[len(output)] = {
                'pos': row_pos,
                'truth_ref': None,
                'truth_alt': None,
                'bcf_ref': row['Ref'],
                'bcf_alt': row['Alt'],
                'snippy_ref': None,
                'snippy_alt': None
            }

    #snippy-uniqe calls
    for i in range(len(snippy_vcf)):
        row = snippy_vcf.iloc[i]
        row_pos = row['Pos']
        if row_pos not in snippy_found:
            output.loc[len(output)] = {
                'pos': row_pos,
                'truth_ref': None,
                'truth_alt': None,
                'bcf_ref': None,
                'bcf_alt': None,
                'snippy_ref': row['Ref'],
                'snippy_alt': row['Alt']
            }


    return output

def without_truth(bcf_vcf, snippy_vcf):
    output = pd.DataFrame(columns=['pos','bcf_ref','bcf_alt','snippy_ref','snippy_alt'])
    window = 10
    snippy_found = [] 

    for i in range(len(bcf_vcf)):
        pos_value = bcf_vcf.at[i, 'Pos']
        bcf_ref = str(bcf_vcf.at[i, 'Ref'])
        bcf_alt = str(bcf_vcf.at[i, 'Alt'])
        snippy_row, snippy_ref, snippy_alt = None, None, None

        #SNP
        if len(bcf_ref) == 1 and len(bcf_alt) == 1:
            snippy_row = snippy_vcf.loc[snippy_vcf['Pos'] == pos_value]
            if not snippy_row.empty:
                snippy_ref = str(snippy_row['Ref'].iloc[0])
                snippy_alt = str(snippy_row['Alt'].iloc[0])
                snippy_found.append(snippy_row['Pos'].iloc[0])

        #insertion
        if len(bcf_ref) < len(bcf_alt):
            min_pos = pos_value - window
            max_pos = pos_value + window
            snippy_row = snippy_vcf.loc[(snippy_vcf['Pos'] >= min_pos) & (snippy_vcf['Pos'] <= max_pos)]
            if not snippy_row.empty:
                snippy_ref = str(snippy_row['Ref'].iloc[0])
                snippy_alt = str(snippy_row['Alt'].iloc[0])
                if len(snippy_alt) > len(snippy_ref):
                    if insertions_match(bcf_ref, bcf_alt, snippy_ref, snippy_alt):
                        snippy_found.append(snippy_row['Pos'].iloc[0])
                    else:
                        snippy_ref = None
                        snippy_alt = None
                else:
                    snippy_ref = None
                    snippy_alt = None

        #deletion
        if len(bcf_ref) > len(bcf_alt):
            min_pos = pos_value - window
            max_pos = pos_value + window
            snippy_row = snippy_vcf.loc[(snippy_vcf['Pos'] >= min_pos) & (snippy_vcf['Pos'] <= max_pos)]
            if not snippy_row.empty:
                snippy_ref = str(snippy_row['Ref'].iloc[0])
                snippy_alt = str(snippy_row['Alt'].iloc[0])
                if len(snippy_ref) > len(snippy_alt):
                    if deletions_match(bcf_ref, bcf_alt, snippy_ref, snippy_alt):
                        snippy_found.append(snippy_row['Pos'].iloc[0])
                    else:
                        snippy_ref = None
                        snippy_alt = None
                else:
                    snippy_ref = None
                    snippy_alt = None

        output.loc[len(output)] = {
            'pos': pos_value,
            'bcf_ref': bcf_ref,
            'bcf_alt': bcf_alt,
            'snippy_ref': snippy_ref,
            'snippy_alt': snippy_alt
        }

    #snippy unique
    for j in range(len(snippy_vcf)):
        row = snippy_vcf.iloc[j]
        row_pos = row['Pos']

        if row_pos not in snippy_found:
            output.loc[len(output)] = {
                'pos': row_pos,
                'bcf_ref': None,
                'bcf_alt': None,
                'snippy_ref': str(row['Ref']),
                'snippy_alt': str(row['Alt'])
            }

    return output

def summary_truth(df, refcol, altcol):
        truth_present = df['truth_ref'].notna() & df['truth_alt'].notna()
        call_present = df[refcol].notna() & df[altcol].notna()

        TP = (truth_present & call_present).sum()
        FP = (~truth_present & call_present).sum()
        FN = (truth_present & ~call_present).sum()

        precision = TP / (TP + FP) if (TP + FP) else 0
        recall    = TP / (TP + FN) if (TP + FN) else 0

        return TP, FP, FN, precision, recall
    
def summary_no_truth(df):
    bcf_present = df['bcf_ref'].notna() & df['bcf_alt'].notna()
    snippy_present = df['snippy_ref'].notna() & df['snippy_alt'].notna()

    unique_bcf = (bcf_present & ~snippy_present).sum()
    unique_snippy = (~bcf_present & snippy_present).sum()
    matches = (bcf_present & snippy_present).sum()

    return matches, unique_bcf, unique_snippy
    
def parse_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bcf", required=True, help="VCF file produced by bcftools.")
    parser.add_argument("-s", "--snippy", required=True, help="VCF file produced by Snippy.")
    parser.add_argument("-t", "--truth", required=False, default=None, help="Optional truth mutations CSV (may or may not exist).")
    parser.add_argument("-o", "--prefix",required=True,help="Prefix for output files (e.g. sample → sample_compare.csv).")

    return parser.parse_args()

def main():
    
    args = parse_args()
    bcf_vcf = args.bcf
    snippy_vcf = args.snippy
    truth_path = args.truth
    prefix = args.prefix
    
    out_df = f"{prefix}_compare.txt"
    out_txt = f"{prefix}_compare_summary.txt"

    truth_exists = truth_path and os.path.exists(truth_path)
    
    f = open(snippy_vcf)
    data = []
    for line in f:
        if line[0] != '#':
            data.append(line.split('\t'))
    snippy_vcf = pd.DataFrame(data, columns=['Chrom', 'Pos', 'ID', 'Ref', 'Alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'file'])

    f = open(bcf_vcf)
    data = []
    for line in f:
        if line[0] != '#':
            data.append(line.split('\t'))
    bcf_vcf = pd.DataFrame(data, columns=['Chrom', 'Pos', 'ID', 'Ref', 'Alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'file'])
    
    snippy_vcf['Pos'] = np.int64(snippy_vcf['Pos'])
    bcf_vcf['Pos'] = np.int64(bcf_vcf['Pos'])

    if truth_exists:
        truth = pd.read_csv(truth_path, index_col=False)
        truth['Pos'] = np.int64(truth['index'] + 1)
        output_df = with_truth(truth, bcf_vcf, snippy_vcf)
        
        bcf_TP, bcf_FP, bcf_FN, bcf_prec, bcf_rec = summary_truth(output_df, 'bcf_ref', 'bcf_alt')
        sn_TP,  sn_FP,  sn_FN,  sn_prec,  sn_rec = summary_truth(output_df, 'snippy_ref', 'snippy_alt')
        
        summary_df = pd.DataFrame([
        {"tool": "bcftools",
        "True Positives": bcf_TP, "False Positives": bcf_FP, "False Negatives": bcf_FN,
        "precision": bcf_prec,
        "recall": bcf_rec,},
        {"tool": "snippy",
        "True Positives": sn_TP, "False Positives": sn_FP, "False Negatives": sn_FN,
        "precision": sn_prec,
        "recall": sn_rec,}])
        
        summary_df.to_csv(out_txt, sep='\t', index=False)
        
    else:
        
        output_df = without_truth(bcf_vcf, snippy_vcf)
        matches, unique_bcf, unique_snippy = summary_no_truth(output_df)
        summary_df = pd.DataFrame([{"matches":matches, "unique bcf":unique_bcf, "unique snippy":unique_snippy}])
        summary_df.to_csv(out_txt, sep='\t', index=False)
        
    
    output_df.to_csv(out_df, sep='\t', index=False)

if __name__ == "__main__":
    main()