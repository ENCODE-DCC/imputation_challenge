#!/usr/bin/env python3
"""Imputation challenge metadata

Author:
    Jin Lee (leepc12@gmail.com)
"""

from logger import log

CELL_NAME = {
    'C02': 'adrenal_gland',
    'C20': 'heart_left_ventricle',
    'C35': 'omental_fat_pad',
    'C45': 'testis',
    'C05': 'BE2C',
    'C06': 'brain_microvascular_endothelial_cell',
    'C07': 'Caco-2',
    'C08': 'cardiac_fibroblast',
    'C11': 'dermis_microvascular_lymphatic_vessel_endothelial_cell',
    'C15': 'G401',
    'C16': 'GM06990',
    'C17': 'H1-hESC',
    'C18': 'H9',
    'C19': 'HAP-1',
    'C21': 'hematopoietic_multipotent_progenitor_cell',
    'C22': 'HL-60',
    'C23': 'IMR-90',
    'C24': 'K562',
    'C27': 'mesenchymal_stem_cell',
    'C28': 'MG63',
    'C30': 'NCI-H460',
    'C32': 'neural_stem_progenitor_cell',
    'C33': 'occipital_lobe',
    'C39': 'SJCRH30',
    'C40': 'SJSA1',
    'C41': 'SK-MEL-5',
    'C42': 'skin_fibroblast',
    'C43': 'skin_of_body',
    'C44': 'T47D',
    'C46': 'trophoblast_cell',
    'C47': 'upper_lobe_of_left_lung',
    'C48': 'urinary_bladder',
    'C49': 'uterus',
    'C51': 'WERI-Rb-1',
    'C12': 'DND-41',
    'C25': 'KMS-11',
    'C31': 'NCI-H929',
    'C34': 'OCI-LY7',
    'C01': 'adipose_tissue',
    'C03': 'adrenal_gland_embryonic',
    'C09': 'CD4-positive_alpha-beta_memory_T_cell',
    'C10': 'chorion',
    'C13': 'endocrine_pancreas',
    'C36': 'peripheral_blood_mononuclear_cell',
    'C37': 'prostate',
    'C38': 'RWPE2',
    'C50': 'vagina',
    'C04': 'amnion',
    'C29': 'myoepithelial_cell_of_mammary_gland',
    'C14': 'ES-I3',
    'C26': 'lower_leg_skin'
}

ASSAY_NAME = {
    'M01': 'ATAC-seq',
    'M02': 'DNase-seq',
    'M03': 'H2AFZ',
    'M04': 'H2AK5ac',
    'M05': 'H2AK9ac',
    'M06': 'H2BK120ac',
    'M07': 'H2BK12ac',
    'M08': 'H2BK15ac',
    'M09': 'H2BK20ac',
    'M10': 'H2BK5ac',
    'M11': 'H3F3A',
    'M12': 'H3K14ac',
    'M13': 'H3K18ac',
    'M14': 'H3K23ac',
    'M15': 'H3K23me2',
    'M16': 'H3K27ac',
    'M17': 'H3K27me3',
    'M18': 'H3K36me3',
    'M19': 'H3K4ac',
    'M20': 'H3K4me1',
    'M21': 'H3K4me2',
    'M22': 'H3K4me3',
    'M23': 'H3K56ac',
    'M24': 'H3K79me1',
    'M25': 'H3K79me2',
    'M26': 'H3K9ac',
    'M27': 'H3K9me1',
    'M28': 'H3K9me2',
    'M29': 'H3K9me3',
    'M30': 'H3T11ph',
    'M31': 'H4K12ac',
    'M32': 'H4K20me1',
    'M33': 'H4K5ac',
    'M34': 'H4K8ac',
    'M35': 'H4K91ac'
}

# cat hg38.chrom.sizes | grep -P "chr[\dX]" | grep -v _
CHRSZ = {
    'chr1': 248956422,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr2': 242193529,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chrX': 156040895
}


def get_cell_name(cell_id):
    if cell_id in CELL_NAME:
        return CELL_NAME[cell_id].replace('_', ' ')
    else:
        return str(cell_id)


def get_assay_name(assay_id):
    if assay_id in ASSAY_NAME:
        return ASSAY_NAME[assay_id].replace('_', ' ')
    else:
        return str(assay_id)


def parse_submission_filename(bw_file):
    """Filename should be CXXMYY.bigwig or CXXMYY.bw
    """
    basename_wo_ext = os.bath.basename(os.path.splitext(bw_file)[0])
    cell = basename_wo_ext[0:3]
    assay = basename_wo_ext[3:6]
    return cell, assay