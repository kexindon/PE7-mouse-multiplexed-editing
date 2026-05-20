"""Run PEGG pegRNA design on all H2M-modelable mouse variants from AACR top-3.

Standalone script for overnight background run.
Resumable: skips chunks whose output CSV already exists.

Usage (from project root, with pegg_env active):
    nohup /opt/miniconda3/envs/pegg_env/bin/python run_pegg_full.py \
        > logs/pegg_full.log 2>&1 &
"""

import os
import sys
import glob
import time
import pandas as pd
from pegg import prime


# === Config ===
CHUNK_SIZE = 5000
OUT_DIR = 'data/pegg_chunks'
INPUT_FORMAT = 'cBioPortal'
PAM = 'NGG'
PBS_LENGTHS = [7]
RTT_LENGTHS = [15]

AACR_MUT = '/Users/kexindong/Documents/GitHub/Database/PublicDatabase/AACR-GENIE/v15.0/data_mutations_extended.txt'
H2M_RESULT = '/Users/kexindong/Documents/GitHub/Output/h2m_database/database-output-final/df_result_cleaned_v3.csv'
PATH_M_REF = '/Users/kexindong/Documents/GitHub/Database/RefGenome/mouse-2023-09-13/GCF_000001635.27_GRCm39_genomic.fna.gz'


def log(msg):
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S")}] {msg}', flush=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs('logs', exist_ok=True)

    # === Load AACR mutations ===
    log('Loading AACR-GENIE mutations...')
    df_mut = pd.read_csv(
        AACR_MUT,
        header=0, sep='\t', comment='#',
        na_values=['Not Applicable', 'NA', 'NULL', 'None', ''],  # '-' kept as valid (indel encoding)
    )
    log(f'  Loaded {len(df_mut):,} mutations')

    # === VAF and top-3 per sample ===
    df_mut['VAF'] = df_mut['t_alt_count'] / (df_mut['t_alt_count'] + df_mut['t_ref_count'])
    df_mut = df_mut[df_mut['VAF'].notna()].reset_index(drop=True)
    log(f'  After VAF filter: {len(df_mut):,}')

    df_mut = (
        df_mut.sort_values(['Tumor_Sample_Barcode', 'VAF'], ascending=[True, False])
              .groupby('Tumor_Sample_Barcode')
              .head(3)
              .reset_index(drop=True)
    )
    log(f'  Top-3 per sample: {len(df_mut):,}')

    # === Unique top-3 mutations for H2M merge ===
    keep_cols = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
                 'Strand', 'Consequence', 'Variant_Classification', 'Variant_Type',
                 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
    df_top3_undup = df_mut[keep_cols].drop_duplicates().reset_index(drop=True)
    log(f'  Unique top-3 mutations: {len(df_top3_undup):,}')

    # Symbol harmonization (match H2M Database naming)
    dict_of_aacr_symbol = {
        'H3F3A':'H3-3A','H3F3B':'H3-3B','WHSC1':'NSD1','GPR124':'ADGRA2','MKL1':'MRTFA',
        'PARK2':'PRKN','MRE11A':'MRE11','MLLT4':'AFDN','CASC5':'KNL1','MEF2BNB-MEF2B':'BORCS8-MEF2B',
        'PAK7':'PAK5','HIST1H1C':'H1-2','HIST1H3B':'H3C1','HIST1H1E':'H1-4','HIST1H2BD':'H2BC5',
        'HIST3H3':'H3-4','HIST1H3D':'H3C2','HIST1H3E':'H3C3','HIST1H3J':'H3C12','HIST1H2BK':'H2BC8',
        'HIST1H3F':'H3C4','HIST1H1B':'H1-5','HIST1H1D':'H1-3','HIST1H2AC':'H2AC1','HIST1H2BJ':'H2BC13',
        'HIST1H4E':'H4C5','FAM46C':'TENT5C','WHSC1L1':'NSD3','SETD8':'KMT5A','LPHN3':'ADGRL3',
        'BAI3':'ADGRB3','SEPT9':'NAPB','BRE':'BABAM2','RFWD2':'COP1','TCEB1':'ELOC','GBA':'GBA1',
        'PVRL4':'NECTIN4','ICK':'GCKR','GNB2L1':'RACK1','MGEA5':'OGA','DIRC2':'SLC49A4',
        'LARGE':'LARGE1','TMEM173B':'STING1','SEPT5':'SEPTIN5'
    }
    df_top3_undup['Hugo_Symbol'] = df_top3_undup['Hugo_Symbol'].replace(dict_of_aacr_symbol)

    # === Load H2M ===
    log('Loading H2M database...')
    db_h2m = pd.read_csv(H2M_RESULT)
    log(f'  Loaded {len(db_h2m):,} H2M rows')

    h2m_lookup = db_h2m[db_h2m['Database'] == 'AACR'].reset_index(drop=True)
    h2m_lookup = h2m_lookup[['gene_name_h', 'chr_h', 'start_h', 'end_h', 'ref_seq_h', 'alt_seq_h',
                              'status', 'ID_human', 'ID_mouse', 'type_h']].copy()
    h2m_lookup['chr_h'] = h2m_lookup['chr_h'].str.replace('^chr', '', regex=True)
    h2m_lookup = h2m_lookup.rename(columns={
        'gene_name_h': 'Hugo_Symbol', 'chr_h': 'Chromosome',
        'start_h': 'Start_Position', 'end_h': 'End_Position',
        'ref_seq_h': 'Reference_Allele', 'alt_seq_h': 'Tumor_Seq_Allele2',
        'status': 'h2m_status', 'type_h': 'Variant_Type',
    })

    # Align dtypes / indel encoding
    df_top3_undup['Chromosome'] = df_top3_undup['Chromosome'].astype(str)
    h2m_lookup['Chromosome'] = h2m_lookup['Chromosome'].astype(str)
    df_top3_undup['Reference_Allele'] = df_top3_undup['Reference_Allele'].fillna('-')
    df_top3_undup['Tumor_Seq_Allele2'] = df_top3_undup['Tumor_Seq_Allele2'].fillna('-')

    # Merge
    df_fused = df_top3_undup.merge(
        h2m_lookup,
        on=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
            'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Type'],
        how='left'
    )
    df_fused['h2m_modelable'] = df_fused['h2m_status'] == True
    n_modelable = df_fused['h2m_modelable'].sum()
    log(f'  H2M-modelable top-3 mutations: {n_modelable:,} / {len(df_fused):,}')

    # === Build PEGG input ===
    list_id_mouse_variant = df_fused[df_fused['h2m_modelable']]['ID_mouse'].dropna().tolist()
    h2m_mouse_for_pegg = db_h2m[db_h2m['ID_mouse'].isin(list_id_mouse_variant)].reset_index(drop=True)
    log(f'  Rows from H2M matching modelable variants: {len(h2m_mouse_for_pegg):,}')

    pegg_input = h2m_mouse_for_pegg.copy()

    # Drop rows where mouse coords are missing (h2m_status=False / no homolog)
    before = len(pegg_input)
    pegg_input = pegg_input.dropna(subset=['start_m', 'end_m', 'chr_m', 'type_m']).reset_index(drop=True)
    log(f'  Rows after dropping missing mouse coords: {len(pegg_input):,} (dropped {before-len(pegg_input):,})')

    pegg_input['Start_Position'] = pegg_input['start_m'].astype(int)
    pegg_input['End_Position'] = pegg_input['end_m'].astype(int)
    pegg_input['Variant_Type'] = pegg_input['type_m']
    pegg_input['Reference_Allele'] = pegg_input['ref_seq_m'].fillna('-')   # INS edge cases have NaN ref
    pegg_input['Tumor_Seq_Allele2'] = pegg_input['alt_seq_m'].fillna('-')  # DEL has '-' alt; safety net
    pegg_input['Chromosome'] = pegg_input['chr_m'].str.replace('chr', '', regex=False)

    # Map MAF variant types -> PEGG vocabulary
    # PEGG only knows: SNP, ONP, DNP, INDEL, INS, DEL
    # MAF uses: SNP, DNP, TNP, ONP -> PEGG treats DNP/TNP/ONP all as multi-bp substitution
    pegg_input['Variant_Type'] = pegg_input['Variant_Type'].replace({
        'TNP': 'ONP',
        'DNP': 'ONP',
    })
    log(f'  Variant type distribution:\n{pegg_input["Variant_Type"].value_counts().to_string()}')

    # Dedup at variant level (drop -V01 suffix)
    pegg_input['ID_mouse_dedup'] = pegg_input['ID_mouse'].str.rsplit('-', n=1).str[0]
    pegg_input = pegg_input.drop_duplicates(subset='ID_mouse_dedup').reset_index(drop=True)
    log(f'  Unique mouse variants for PEGG: {len(pegg_input):,}')

    # === Load mouse genome ===
    log('Loading mouse reference genome (this takes a few minutes)...')
    t_gen = time.time()
    chrom_dict, _ = prime.genome_loader(PATH_M_REF)
    log(f'  Genome loaded in {time.time()-t_gen:.0f}s')

    # === Chunked PEGG run ===
    n_chunks = (len(pegg_input) + CHUNK_SIZE - 1) // CHUNK_SIZE
    log(f'\nRunning {n_chunks} chunks of size {CHUNK_SIZE:,}')

    overall_start = time.time()
    for i in range(n_chunks):
        out_path = f'{OUT_DIR}/pegRNAs_chunk_{i:03d}.csv'
        if os.path.exists(out_path):
            log(f'[{i+1:3d}/{n_chunks}] skipped (already exists)')
            continue

        t0 = time.time()
        sub = pegg_input.iloc[i*CHUNK_SIZE:(i+1)*CHUNK_SIZE].copy()
        try:
            pegRNAs_chunk = prime.run(
                sub,
                INPUT_FORMAT,
                chrom_dict=chrom_dict,
                PAM=PAM,
                PBS_lengths=PBS_LENGTHS,
                RTT_lengths=RTT_LENGTHS,
            )
            pegRNAs_chunk.to_csv(out_path, index=False)
            elapsed = time.time() - t0
            log(f'[{i+1:3d}/{n_chunks}] {len(sub):,} muts -> {len(pegRNAs_chunk):,} pegRNAs '
                f'in {elapsed:.0f}s ({len(sub)/elapsed:.2f} muts/s)')
        except Exception as e:
            log(f'[{i+1:3d}/{n_chunks}] FAILED: {e}')
            continue

    log(f'\nAll chunks done in {(time.time()-overall_start)/3600:.2f} hours')

    # === Combine ===
    chunk_files = sorted(glob.glob(f'{OUT_DIR}/pegRNAs_chunk_*.csv'))
    if not chunk_files:
        log('No chunk files found; nothing to combine.')
        return

    pegRNAs_full = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)
    log(f'Combined pegRNAs: {len(pegRNAs_full):,}')

    pegRNAs_filt = pegRNAs_full[pegRNAs_full['contains_polyT_terminator'] == False]
    pegRNAs_filt = pegRNAs_filt[
        pegRNAs_filt['sensor_error'] != 'sensor too small; increase sensor size in parameters; or decrease before_proto_context'
    ]
    pegRNAs_filt.to_csv('data/pegRNAs_full_filtered.csv', index=False)
    log(f'Filtered pegRNAs: {len(pegRNAs_filt):,}')

    coord_cols = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
    n_amenable = pegRNAs_filt[coord_cols].drop_duplicates().shape[0]
    log(f'Unique PE7-amenable mouse variants: {n_amenable:,}')
    log('DONE')


if __name__ == '__main__':
    main()
