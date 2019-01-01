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

import argparse
import gzip


####################################################################################################

def parse_pfam2go(pfam2go='pfam2go.txt'):
    """
    :param pfam2go: full path to a file containing Pfam->GO annotations
    :return: Pfam ID -> set(binding it is associated with)
    """

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
# MAIN
####################################################################################################

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Highlight somatically perturbed mechanisms in putative cancer'
                                                 'driver genes.')
    parser.add_argument('--pertinint_results', type=str, help='Path to output file from PertInInt', default=None)
    parser.add_argument('--pfam2go_file', type=str, help='Path to Pfam2Go file (will be downloaded if not present)',
                        default='pfam2go.txt')
    parser.add_argument('--out_file', type=str, help='Results output file', default=None)
    args = parser.parse_args()
