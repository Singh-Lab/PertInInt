#!/usr/bin/python

"""
PertInInt:
Analytical approach to rapidly uncover proteins with significant enrichments of somatic Perturbations In Interaction
and other functional sites.

Please cite:
Kobren, S.N., Chazelle, B. and Singh, M. (2019). "An integrative approach to identify preferentially altered
interactions in human cancers." Manuscript in preparation.

GitHub page:
https://github.com/Singh-Lab/PertInInt

Please contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) or Mona Singh (mona@cs.princeton.edu) with questions.
"""

import os
import sys
import argparse
import gzip
from subprocess import call


####################################################################################################
# CONSTANTS
####################################################################################################

common_false_positives = {
    'ABCA1', 'ABCA12', 'ABCB1', 'ABCC9', 'ABCD2', 'ABCG2', 'ABLIM3', 'ACACA', 'ACAN', 'ACSM5',
    'ACTN2', 'ADAM19', 'ADAM23', 'ADAMTS12', 'ADAMTS16', 'ADAMTS18', 'ADAMTS20', 'ADAMTS8',
    'ADAMTSL1', 'ADAMTSL3', 'ADCY1', 'ADCY2', 'ADCY5', 'AGAP10', 'AKAP6', 'AMPH', 'ANK1',
    'ANKRD11', 'ANKS1B', 'ANO3', 'ANTXR1', 'AOC3', 'ARAP2', 'ARAP3', 'ARFGEF1', 'ARHGAP5',
    'ARMC3', 'ASH1L', 'ASTN2', 'ASXL3', 'ATP11C', 'ATP1A2', 'ATP8A2', 'ATP8B4', 'ATXN1', 'AUTS2',
    'BAI3', 'BCHE', 'BCL11A', 'BCL11B', 'BCL9', 'BCLAF1', 'BICC1', 'BIRC6', 'BMPER', 'BNC2',
    'BRINP1', 'BRINP3', 'BRWD3', 'BTBD11', 'BTK', 'BTRC', 'C1S', 'CACHD1', 'CACNA2D1', 'CACNA2D3',
    'CACNB2', 'CADPS', 'CADPS2', 'CALCR', 'CAMTA1', 'CAPN6', 'CARD8', 'CASZ1', 'CCDC39', 'CDC27',
    'CDC73', 'CDH11', 'CDH12', 'CDH2', 'CDH8', 'CEP170', 'CHD7', 'CHL1', 'CHRM3', 'CLCA2',
    'CLIP1', 'CLSTN2', 'CNTN1', 'CNTN3', 'CNTNAP2', 'CNTNAP3B', 'COL12A1', 'COL14A1', 'COL1A2',
    'COL21A1', 'COL25A1', 'COL4A5', 'COLGALT2', 'CPED1', 'CPNE4', 'CPS1', 'CPXM2', 'CR2', 'CSMD1',
    'CSMD2', 'CTCFL', 'CTNNA2', 'CTNND2', 'CUBN', 'CYP26B1', 'DCC', 'DCLK1', 'DDR2', 'DDX26B',
    'DGKI', 'DLG2', 'DLGAP3', 'DNAH11', 'DNAH3', 'DNAH5', 'DNAH7', 'DNAH8', 'DOCK2', 'DOCK4',
    'DPP10', 'DPP6', 'DPYS', 'DSCAML1', 'DTNA', 'DTX1', 'DYNC1I1', 'DYNC2H1', 'DYSF', 'DZIP3',
    'EBF1', 'EBF3', 'EGFLAM', 'ELMO1', 'EP400', 'EPHA5', 'EPHA6', 'EPHA7', 'EPHB1', 'ESRRG',
    'EXOC2', 'EXOC4', 'EYA1', 'EYA4', 'F13A1', 'FAM46A', 'FAM83B', 'FANCD2', 'FAT2', 'FBN2',
    'FBXL7', 'FGGY', 'FHOD3', 'FIGN', 'FLG', 'FLNC', 'FLT1', 'FMN2', 'FMO1', 'FOLH1', 'FOXC2',
    'FOXP2', 'FOXQ1', 'FREM2', 'FYB', 'FZD10', 'FZD8', 'GLCCI1', 'GLDC', 'GLG1', 'GLI1', 'GLI2',
    'GLI3', 'GPC6', 'GPR125', 'GPR158', 'GPRASP2', 'GRIA3', 'GRID1', 'GRIK2', 'GRIK3', 'GRIK4',
    'GRIK5', 'GRIN2A', 'GRIN3A', 'GRM5', 'GRM6', 'GRM8', 'GUCY1A2', 'HAS2', 'HCN1', 'HDX',
    'HECTD2', 'HECTD4', 'HERC2', 'HHIPL2', 'HMCN1', 'HNF1A', 'HRNR', 'HSPA1L', 'HSPA6', 'IGSF9B',
    'IKZF1', 'IL32', 'IL7R', 'INPP4B', 'INPP5D', 'IQGAP2', 'ITGA10', 'ITGA4', 'ITGA8', 'ITGAV',
    'ITIH5', 'ITPR2', 'ITSN1', 'KCND2', 'KCNG1', 'KCNH4', 'KCNT2', 'KIAA1109', 'KIAA1211',
    'KIAA1324L', 'KIAA2022', 'KIF21B', 'KIF4B', 'KIF5A', 'KIRREL', 'KLHL4', 'KLHL5', 'LAMA2',
    'LAMA4', 'LAMB1', 'LGI2', 'LGR5', 'LINGO1', 'LPHN2', 'LPHN3', 'LRFN5', 'LRP12', 'LRP2',
    'LRP4', 'LRP6', 'LRRC4B', 'LRRIQ1', 'MADD', 'MAGI2', 'MAP1B', 'MARK1', 'MED12L', 'MEF2A',
    'MEGF10', 'MGAM', 'MGAT3', 'MGAT5B', 'MLH1', 'MLLT10', 'MME', 'MMP16', 'MMP2', 'MUC16',
    'MUC6', 'MXRA5', 'MYB', 'MYCBP2', 'MYH3', 'MYO3B', 'MYOM2', 'NAALAD2', 'NALCN', 'NBAS',
    'NBEA', 'NCAM1', 'NCAM2', 'NCKAP5', 'NCOR2', 'NEB', 'NELL2', 'NLGN1', 'NLRP3', 'NOTCH2',
    'NOVA1', 'NPAS3', 'NPR3', 'NR4A2', 'NRCAM', 'NRG3', 'NRK', 'NRXN1', 'NRXN3', 'NTNG1', 'NTRK2',
    'NYAP2', 'OCA2', 'OGDHL', 'OGT', 'OPHN1', 'OTOF', 'OTOGL', 'OTUD7A', 'OVGP1', 'PABPC3',
    'PAIP1', 'PAK3', 'PAMR1', 'PAPPA', 'PARD3B', 'PCDH10', 'PCDH11X', 'PCDH17', 'PCDH18',
    'PCDH19', 'PCDH7', 'PCDH9', 'PCDHB16', 'PCDHB17', 'PCDHB8', 'PCDHGA2', 'PCDHGA3', 'PCDHGA6',
    'PCLO', 'PCMTD1', 'PCNX', 'PCSK2', 'PCSK5', 'PDE10A', 'PDE11A', 'PDE1A', 'PDE1C', 'PDE3A',
    'PDZD2', 'PDZRN4', 'PEG3', 'PHACTR4', 'PHEX', 'PIAS3', 'PIK3C2B', 'PKHD1L1', 'PLA2G4A',
    'PLCB1', 'PLCB4', 'PLEC', 'PLEKHA7', 'PLXNA4', 'PODXL', 'POLR3B', 'POSTN', 'POT1', 'POTEE',
    'POTEF', 'PPARGC1A', 'PPFIA2', 'PPP1R10', 'PPP1R16B', 'PPP1R9A', 'PPP3CA', 'PREX1', 'PREX2',
    'PRKD1', 'PRKG1', 'PROS1', 'PROX1', 'PTCHD2', 'PTCHD4', 'PTPRB', 'PTPRC', 'PTPRZ1', 'PXDN',
    'RAG1', 'RALGAPB', 'RALY', 'RANBP6', 'RARG', 'RASAL2', 'RGPD3', 'RGPD4', 'RIMBP2', 'RIMS1',
    'RIMS2', 'ROBO1', 'ROBO2', 'ROR2', 'RUNX1T1', 'RYR2', 'RYR3', 'SACS', 'SALL2', 'SATB2',
    'SCN1A', 'SCN2A', 'SCN5A', 'SCN9A', 'SEC24C', 'SEMA3A', 'SEMA3D', 'SEMA5A', 'SEMA6A',
    'SEMA6D', 'SFMBT2', 'SIPA1L1', 'SIRPA', 'SKIDA1', 'SLC12A5', 'SLC26A7', 'SLC44A5', 'SLC4A10',
    'SLC4A4', 'SLC8A1', 'SLC8A3', 'SLC9A2', 'SLCO5A1', 'SLIT3', 'SLITRK6', 'SPEF2', 'SPEG',
    'SPEN', 'SPTA1', 'SRCAP', 'SRGAP3', 'SSH2', 'ST18', 'ST6GAL2', 'STK19', 'SULF1', 'SUPT16H',
    'SVEP1', 'SYNE1', 'TAF1L', 'TAF3', 'TARS', 'TARS2', 'TBC1D22A', 'TBX5', 'TCF20', 'TCF4',
    'TECTA', 'TENM1', 'TG', 'THBS2', 'THOC2', 'THSD7A', 'THSD7B', 'TMEM132B', 'TMTC2', 'TOX2',
    'TP63', 'TRIO', 'TRPA1', 'TRPC4', 'TRPC6', 'TRPM8', 'TRPS1', 'TRRAP', 'TXNIP', 'UBC', 'UBE3A',
    'UHRF1BP1', 'UHRF1BP1L', 'UNC45B', 'USH2A', 'VANGL2', 'VAV3', 'VPS13D', 'VPS52', 'WAC',
    'WBSCR17', 'WDFY3', 'WDR70', 'WSCD2', 'ZBTB20', 'ZC3H12B', 'ZEB2', 'ZFHX4', 'ZFPM2', 'ZNF148',
    'ZNF292', 'ZNF423', 'ZNF521', 'ZNF536', 'ZNF626', 'ZNF706', 'ZNF800', 'ZNF831', 'ZNF91'
}


####################################################################################################
# BINDING ANNOTATION MAPPINGS
####################################################################################################

def download_pfam2go(pfam2go='pfam2go.txt'):
    """
    :param pfam2go: full path to the Pfam2Go file
    :return: full path to the Pfam2Go file
    """

    # text file of Pfam IDs -> GO annotations
    if not os.path.isfile(pfam2go):
        call(['wget',
              'http://www.geneontology.org/external2go/pfam2go',
              '-O',
              pfam2go])

    if not os.path.isfile(pfam2go):
        sys.stderr.write('Failed to download Pfam2Go annotations from http://www.geneontology.org/external2go/' +
                         'pfam2go... exiting\n')
        sys.exit(1)

    return pfam2go


####################################################################################################

def download_ligand_mapping(ligand_grps='interacdome_ligandgrps.txt'):
    """
    :param ligand_grps: full path to the ligand mmCIF ID -> ligand group mapping file
    :return: full path to the ligand mmCIF ID -> ligand group mapping file
    """

    # tab-delimited text file of ligand type -> mmCIF ID
    if not os.path.isfile(ligand_grps):
        call(['wget',
              'https://raw.githubusercontent.com/Singh-Lab/InteracDome/master/downloaded_data/ligand_groups.txt',
              '-O',
              ligand_grps])

    if not os.path.isfile(ligand_grps):
        sys.stderr.write('Failed to download mmCIF -> ligand group mapping from https://raw.githubusercontent.com/' +
                         'Singh-Lab/InteracDome/master/downloaded_data/ligand_groups.txt... exiting\n')
        sys.exit(1)

    return ligand_grps


####################################################################################################

def parse_pfam2go(pfam2go='pfam2go.txt'):
    """
    :param pfam2go: full path to a file containing Pfam->GO annotations
    :return: Pfam ID -> set(binding types it is associated with)
    """

    if not pfam2go or not os.path.isfile(pfam2go):
        pfam2go = download_pfam2go(pfam2go)

    pfam_to_binding = {}

    precedence_order = ['domain_signaling', 'domain_rna', 'domain_peptide', 'domain_dna', 'domain_metabolite',
                        'domain_ion']

    pfam_handle = gzip.open(pfam2go) if pfam2go.endswith('gz') else open(pfam2go)
    for pfam_line in pfam_handle:
        if pfam_line.startswith('!'):
            continue
        pfam_id = pfam_line[5:12] + '_' + pfam_line.split()[1]  # e.g., "PF00096_zf-C2H2"
        pfam_desc = pfam_line[pfam_line.find('>') + 1:pfam_line.rfind(';')].replace('GO:', '').strip()

        if 'signaling' in pfam_desc or 'pathway' in pfam_desc:
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_signaling') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_signaling'

        elif 'DNA' in pfam_desc and 'binding':
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_dna') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_dna'

        elif 'RNA binding' in pfam_desc:
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_rna') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_rna'

        elif 'protein' in pfam_desc and 'binding' in pfam_desc:
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_peptide') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_peptide'

        elif 'ion binding' in pfam_desc:
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_ion') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_ion'

        elif 'binding' in pfam_desc and True not in [lt in pfam_desc for lt in ['DNA', 'RNA', 'protein', 'ion binding',
                                                                                'nucleic acid binding']]:
            if pfam_id not in pfam_to_binding or \
               precedence_order.index('domain_metabolite') < precedence_order.index(pfam_to_binding[pfam_id]):
                pfam_to_binding[pfam_id] = 'domain_metabolite'

    pfam_handle.close()

    return pfam_to_binding


####################################################################################################

def parse_ligandgrps(ligand_grps='interacdome_ligandgrps.txt'):
    """
    :param ligand_grps: full path to a tab-delimited file containing ligand group to mmCIF ID(s)
    :return: mmCIF -> ligand group type
    """

    if not ligand_grps or not os.path.isfile(ligand_grps):
        ligand_grps = download_ligand_mapping(ligand_grps)

    map_ligands = {'DNA_': 'dna', 'DNABASE_': 'dna', 'DNABACKBONE_': 'dna',
                   'RNA_': 'rna', 'RNABASE_': 'rna', 'RNABACKBONE_': 'rna',
                   'DRUGLIKE_': 'sm', 'METABOLITE_': 'sm', 'SM_': 'sm',  # 'NUCACID_': 'dna',
                   'PEPTIDE_': 'peptide'}

    ligand_handle = gzip.open(ligand_grps) if ligand_grps.endswith('gz') else open(ligand_grps)
    for lig_line in ligand_handle:
        if lig_line.startswith('#'):
            continue

        ligand_group, original_id = lig_line.strip().split('\t')[:2]
        if ligand_group == 'ION_' and original_id not in map_ligands:
            map_ligands[original_id] = 'ion'
        elif ligand_group in ['METABOLITE_', 'DRUGLIKE_', 'SM_'] and original_id not in map_ligands:
            map_ligands[original_id] = 'sm'
    ligand_handle.close()

    return map_ligands


####################################################################################################
# PARSE PERTININT RESULTS
####################################################################################################

def parse_tracknames(track_names, map_ligands, pfam_to_binding):
    """
    :param track_names: ';'-delimited string corresponding to track names and positive zscores
    :param map_ligands: dictionary of mmCIF ID -> ligand grouping
    :param pfam_to_binding: dictionary of Pfam ID -> interaction type
    :return: dictionary of mechanism type -> zscore
    """

    mech_scores = {}
    try:
        track_to_zscore = {tz.split('|')[0]: float(tz.split('|')[1]) for tz in track_names.split(';')}
    except IndexError:
        print track_names
        return mech_scores

    for track_name, zscore in track_to_zscore.items():

        classifications = set()

        for bd_track in track_name.split(','):
            current_classification = ''

            if 'ConCavity' in bd_track:
                current_classification = 'pocket'

            elif 'JSD_conservation' in bd_track:
                current_classification = 'conservation'

            elif 'WholeGene_NatVar' in bd_track:
                current_classification = 'mutfreq'

            elif bd_track.startswith('PF') and bd_track.endswith(':complete'):
                current_classification = pfam_to_binding.get(bd_track[:bd_track.find(':')], 'domain')

            elif bd_track.startswith('Homology_binding') or bd_track.startswith('PF'):
                ligand = bd_track.split(':')[-1]
                if ligand.endswith('_') and ligand not in map_ligands:
                    continue
                current_classification = 'interaction_' + map_ligands.get(ligand,
                                                                          'sm')  # ligand.replace('_','').lower()

            classifications.add(current_classification)

        for track_class in classifications:
            if track_class not in mech_scores:
                mech_scores[track_class] = zscore
            else:
                mech_scores[track_class] = max(zscore, mech_scores[track_class])

    return mech_scores  # mechanism type -> zscore


####################################################################################################

def parse_ordered_genes(results_file, output_file, pfam2go='pfam2go.txt', ligand_grps='interacdome_ligandgrps.txt'):
    """
    :param results_file: full path to a tab-delimited PertInInt results file
    :param output_file: full path to a file to write tab-delimited mechanism output to
    :param pfam2go: full path to a Pfam -> GO annotation mapping file
    :param ligand_grps: full path to a mmCIF -> ligand type mapping file
    :return: full path to successfully written output file
    """

    # (1) get the mapping from Pfam domain ID and ligand type to interaction types:
    map_ligands = parse_ligandgrps(ligand_grps)  # ligand type -> binding
    pfam_to_binding = parse_pfam2go(pfam2go)  # domain -> binding

    # (2) track types corresponding to "mechanisms"
    mechanism_types = ['interaction_dna', 'interaction_rna', 'interaction_peptide',
                       'interaction_ion', 'interaction_metabolite', 'interaction_sm',
                       'domain_dna', 'domain_rna', 'domain_peptide', 'domain_metabolite',
                       'domain_signaling', 'domain', 'conservation', 'mutfreq']

    # (3) start output file
    out_handle = gzip.open(output_file, 'wt') if output_file.endswith('gz') else open(output_file, 'w')
    out_handle.write('\n'.join([
        '# Somatically perturbed functional mechanisms uncovered by',
        '# PertInInt, v0: https://github.com/Singh-Lab/PertInInt',
        '# Kobren, S.N., Chazelle, B. and Singh, M. (2019) "An integrative approach to ' +
        'identify preferentially altered interactions in human cancers." ' +
        'Manuscript in preparation.',
        '# PertInInt results file: ' + results_file,
        '# Pfam2Go annotations: ' + pfam2go,
        '# mmCIF ligand ID annotations: ' + ligand_grps,
        '\t'.join(['cancer_status', 'gene_name'] + mechanism_types)
    ])+'\n')

    # (4) parse input file and write to output file
    in_handle = gzip.open(results_file) if results_file.endswith('gz') else open(results_file)
    header = None
    for result in in_handle:
        if result.startswith('#'):
            continue
        elif not header:
            header = result[:-1].split('\t')
            continue

        v = result[:-1].split('\t')
        cancer_status = v[header.index('cancer_status')]
        gene_name = v[header.index('gene_name')]
        track_names = v[header.index('zscores_per_track')]

        track_types = parse_tracknames(track_names, map_ligands, pfam_to_binding)

        zscores = {}
        for tclass, zscore in track_types.items():
            if tclass not in mechanism_types:
                continue
            zscores[tclass] = zscore

        out_handle.write('\t'.join([cancer_status, gene_name] +
                                   [zscores.get(tclass, '--') for tclass in mechanism_types])+'\n')
    in_handle.close()
    out_handle.close()

    return output_file


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Highlight somatically perturbed mechanisms in putative cancer'
                                                 'driver genes.')
    parser.add_argument('--pertinint_results', type=str, help='Path to output file from PertInInt', default=None)
    parser.add_argument('--out_file', type=str, help='Results output file', default=None)
    parser.add_argument('--pfam2go_file', type=str, help='Path to Pfam2Go file (will be downloaded if not present)',
                        default='pfam2go.txt')
    parser.add_argument('--ligand_groups_file', type=str, help='Path to ligand ID -> ligand group mapping file ' +
                        '(will be downloaded if not present)', default='interacdome_ligandgrps.txt')
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------
    # (i) confirm that we can open PertInInt results file
    if not args.pertinint_results or not os.path.isfile(args.pertinint_results):
        sys.stderr.write('Could not open PertInInt results file: '+str(args.pertinint_results)+'\n' +
                         'Usage: python '+sys.argv[0]+' --pertinint_results <results_file> --out_file <output_file>\n')
        sys.exit(1)

    # (ii) confirm that we can write to mechanisms output file
    if not args.out_file:
        sys.stderr.write('Could not write to specified output file: ' + str(args.out_file) + '\n' +
                         'Usage: python '+sys.argv[0]+' --pertinint_results <results_file> --out_file <output_file>\n')
        sys.exit(1)

    for subdir in ['/'.join(args.out_file.split('/')[:ind]) for ind in xrange(2, args.out_file.count('/')+1)]:
        if not os.path.isdir(subdir):
            if call(['mkdir', subdir]):  # any code returned other than "0"
                sys.stderr.write('Could not write to '+args.out_file+'. Exiting\n')
                sys.exit(1)

    # ------------------------------------------------------------------------------------------------
    parse_ordered_genes(args.pertinint_results,
                        args.out_file,
                        args.pfam2go_file,
                        args.ligand_groups_file)
    sys.stderr.write('Wrote tab-delimited mechanisms to: '+args.out_file+'\n')
