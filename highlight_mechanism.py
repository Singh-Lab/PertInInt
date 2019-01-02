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
        cancer_status = v[header.index('cancer_status')].strip()
        gene_name = v[header.index('gene_name')].strip()
        track_names = v[header.index('zscores_per_track')].strip()

        if len(track_names) < 1:
            continue

        track_types = parse_tracknames(track_names, map_ligands, pfam_to_binding)

        zscores = {}
        for tclass, zscore in track_types.items():
            if tclass not in mechanism_types:
                continue
            zscores[tclass] = zscore

        out_handle.write('\t'.join([cancer_status, gene_name] +
                                   [str(zscores.get(tclass, '--')) for tclass in mechanism_types])+'\n')
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
