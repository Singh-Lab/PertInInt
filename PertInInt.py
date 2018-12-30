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
import time
import gzip
import argparse
import numpy as np
from math import sqrt
from subprocess import call


####################################################################################################
# CONSTANTS
####################################################################################################

RESTRICTED_DOMS = None  # if there are particular (domain, ligand) pairs you aren't interested in, specify them here


####################################################################################################

def reformat_time(run_time):
    """
    :param run_time: total time elapsed in seconds
    :return: string corresponding to a properly formatted (days, hours, minutes, seconds) time
    """

    m, s = divmod(run_time, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    return ':'.join(map(lambda x: str(int(x)).zfill(2), [d, h, m, s]))


####################################################################################################
# FIND, FORMAT AND CHECK TRACKS
####################################################################################################

def track_weights_list(trackweight_path):
    """
    :param trackweight_path: full path to directories containing weight vectors
    :return: protein ID -> full path to track weight file list
    """

    protid_to_track = {}  # Ensembl protein ID -> full path to track file
    protid_to_geneid = {}  # set of modelable Ensembl gene IDs
    suffix = '.trackweights.tsv'  # properly formatted track weight files will end with this suffix

    for chrom_id in [d for d in os.listdir(trackweight_path) if os.path.isdir(trackweight_path+d)]:
        for gene_id in [d for d in os.listdir(trackweight_path+chrom_id)
                        if os.path.isdir(trackweight_path+chrom_id+'/'+d)]:
            for track_file in [d for d in os.listdir(trackweight_path+chrom_id+'/'+gene_id) if d.endswith(suffix)]:

                prot_id = track_file[:track_file.find(suffix)]
                protid_to_track[prot_id] = trackweight_path + chrom_id + '/' + gene_id + '/' + track_file
                protid_to_geneid[prot_id] = gene_id

    return protid_to_track, protid_to_geneid


####################################################################################################

def track_name_to_classification(track_name):
    """
    :param track_name: string corresponding to a track name (from a weight vector file)
    :return: the classification of that track
    """

    # highly correlated tracks were merged into the same track in the preprocess step. It is HIGHLY unlikely
    # that tracks from different categories were merged, but in this rare case, we pick the more "mechanistically
    # informative" classification:
    for bd_track in track_name.split(','):
        if bd_track.startswith('Homology_binding') or \
            (bd_track.startswith('PF') and not bd_track.endswith(':complete') and
             ((not RESTRICTED_DOMS and bd_track.startswith('PF') and not bd_track.endswith(':complete')) or
              (RESTRICTED_DOMS and (bd_track.split(':')[0], bd_track.split(':')[-1]) not in RESTRICTED_DOMS))):
            return 'interaction'
        elif bd_track.startswith('PF') and bd_track.endswith(':complete'):
            return 'domain'
        elif bd_track.startswith('JSD_conservation'):
            return 'conservation'
        elif bd_track.startswith('ExAC_allelefreq'):
            return 'conservation'
        elif bd_track.startswith('ConCavity'):
            return 'interaction'
    return 'none'


####################################################################################################

def edit_track_name(track_name):
    """
    :param track_name: string corresponding to a track name (as in weight vector files)
    :return: reformatted track_name (for reporting mechanism)
    """

    edited_tracks = []
    for bd_track in track_name.split(','):
        if bd_track.startswith('PF') and \
           not bd_track.endswith(':complete') and \
           RESTRICTED_DOMS and (bd_track.split(':')[0], bd_track.split(':')[-1]) in RESTRICTED_DOMS:
            continue
        elif bd_track.startswith('PF') and bd_track.split(':')[1] != 'aggregate':
            edited_tracks.append(bd_track[:bd_track.find(':')] + ':X:' + bd_track[bd_track.rfind(':') + 1:])
        else:
            edited_tracks.append(bd_track)

    return ','.join(edited_tracks)


####################################################################################################

def check_restrictions(track_name, restriction='none', aggregate_names=None):
    """
    :param track_name: string from a weight vector file describing the current weight track
    :param restriction: string indicating a subset of tracks to consider
    :param aggregate_names: names of domains that appear repetitively (in case we want to consider only
                            the aggregate case rather than individual domains)
    :return: boolean whether the track_name passes the specified restriction or not
    """

    # Don't pass any interaction domain or domain track that has a singular member of a repeat domain family
    if aggregate_names:
        for bd_track in track_name.split(','):
            if bd_track.startswith('PF') and bd_track[:bd_track.find(':')] in aggregate_names and \
               bd_track.split(':')[1] != 'aggregate':
                return False  # don't consider individual domains in repetitive families

    track_classification = track_name_to_classification(track_name)
    if track_classification not in ['interaction', 'domain', 'conservation', 'wholegene']:
        return False

    if restriction in ['interaction', 'interwholegene']:
        if track_classification != 'interaction':
            return False
    elif restriction in ['domain', 'domwholegene']:
        if track_classification != 'domain':
            return False
    elif restriction in ['conservation', 'conswholegene']:
        if track_classification != 'conservation':
            return False
    elif restriction in ['nointeraction', 'domcons']:
        if track_classification == 'interaction':
            return False
    elif restriction in ['nodomain', 'intercons']:
        if track_classification == 'domain':
            return False
    elif restriction in ['noconservation', 'interdom']:
        if track_classification == 'conservation':
            return False

    return True  # if we passed, or if our restriction is 'none' or 'nowholegene'


####################################################################################################

def aggregate_domain_tracks(weightfile, repeat_domain_limit=40):
    """
    :param weightfile: full path to a weight vector file
    :param repeat_domain_limit: minimum number of times a domain must be repeated to replace ALL
                                individual instances with the aggregate track only
    :return: set of domains that we found more than X times
    """

    domain_tracks = []
    weightfile_handle = gzip.open(weightfile) if weightfile.endswith('gz') else open(weightfile)
    for wline in weightfile_handle:
        if wline.startswith('#'):
            continue

        track_name = wline[:-1].split('\t')[2]
        domain_tracks.extend([bd_track.split(':')[0] for bd_track in track_name.split(',')
                              if bd_track.startswith('PF') and
                              bd_track.endswith(':complete') and
                              not bd_track.endswith(':aggregate:complete')])  # e.g., PF00096_zf-C2H2:1-23:complete
    weightfile_handle.close()

    return [dname for dname in set(domain_tracks) if domain_tracks.count(dname) >= repeat_domain_limit]


####################################################################################################
# PROCESS MUTATIONS
####################################################################################################

def gene_name_mapping(mapping_file):
    """
    :param mapping_file: full path to a tab-delimited file with Ensembl gene ID in the first column and
                         a comma-delimited list of all gene names and synonyms in the third column
    :return: dictionary of gene name -> ensembl gene ID
    """

    name_mapping = {}  # gene name -> ensembl gene ID

    mapping_handle = gzip.open(mapping_file) if mapping_file.endswith('gz') else open(mapping_file)
    for mapping_line in mapping_handle:
        if mapping_line.startswith('#'):
            continue
        ensembl_id, _, alt_ids = mapping_line[:-1].split('\t')[:3]
        for alt_id in alt_ids.split(','):
            name_mapping[alt_id] = ensembl_id
    mapping_handle.close()

    return name_mapping


####################################################################################################

def process_mutations_from_maf(maf_file, modelable_genes, modelable_prots, mapping_file, expression_file):
    """
    :param maf_file: full path to a .maf formatted file with somatic mutations
    :param modelable_genes: dictionary of Ensembl gene IDs -> modelable protein IDs
    :param modelable_prots: set of Ensembl protein IDs that we want to collect mutation information for
    :param mapping_file: full path to a tab-delimited file with Ensembl gene ID in the first column and
                         a comma-delimited list of all gene names and synonyms in the third column
    :param expression_file: full path to a tab-delimited file containing Ensembl gene IDs and the set of tumor
                            samples in which that gene is expressed
    :return: (1) dictionary prot_id -> [(0-index position of missense mutation, mutational value), ...]
             (2) dictionary prot_id -> total nonsynonymous mutational burden (allelic fraction)
             (3) total nonsynonymous mutational burden across all genes
             (4) total nonsynonymous mutation count across all genes
    """

    # ------------------------------------------------------------------------------------------------
    # (1) get mapping from ensembl gene ID -> set of tumor samples in which this gene is expressed
    expression_by_gene = None
    if expression_file and mapping_file:
        expression_by_gene = {}
        exp_handle = gzip.open(expression_file) if expression_file.endswith('gz') else open(expression_file)
        for exp_line in exp_handle:
            if exp_line.startswith('#'):
                continue
            gene_name = exp_line.split('\t')[0].split(',')[0].split('.')[0]  # ensembl gene ID (no version)
            expression_by_gene[gene_name] = set(exp_line.split('\t')[1].split(','))
        exp_handle.close()

    # (2) get mapping from gene "name" -> ensembl gene ID
    name_to_ensembl = gene_name_mapping(mapping_file)

    # ------------------------------------------------------------------------------------------------
    # (3) define empty variables (to return)
    mutation_values = {prot_id: 0. for prot_id in modelable_prots}  # prot -> total mutational burden
    mutation_locations = {prot_id: [] for prot_id in modelable_prots}  # prot -> [(loc, type, value), ...]
    total_mutational_value = 0.
    total_mutations = 0

    # ------------------------------------------------------------------------------------------------
    # (3) start reading mutation .maf file
    header = None
    mut_handle = gzip.open(maf_file) if maf_file.endswith('gz') else open(maf_file)
    for mutline in mut_handle:
        if mutline.startswith('#'):
            continue
        if not header:
            header = mutline.lower()[:-1].split('\t')
            continue

        # get contents for this line:
        v = mutline[:-1].split('\t')
        prot_id = v[header.index('ensp')]
        gene_name = v[header.index('hugo_symbol')]
        ensembl_id = name_to_ensembl.get(gene_name, None)
        mut_type = v[header.index('variant_classification')].replace('_Mutation', '')
        sample_id = '-'.join(v[header.index('tumor_sample_barcode')].split('-')[:4])
        mut_val = float(v[header.index('t_alt_count')]) / float(v[header.index('t_depth')])

        # make sure that this mutation occurred within a protein:
        try:
            aachange = mutline[:-1].split('\t')[header.index('hgvsp_short')][2:]  # e.g., p.L414L
            mut_pos = int(''.join([i for i in list(aachange) if i in map(str, range(10))])) - 1
        except ValueError:
            continue

        # make sure this mutation is occurring in a gene that is expressed
        if expression_by_gene and sample_id not in expression_by_gene.get(ensembl_id, []):
            continue

        # keep track of all nonsynonymous mutations in modelable genes/proteins
        if mut_type in ['Missense', 'Nonsense'] and ensembl_id in modelable_genes:
            total_mutations += 1
            total_mutational_value += mut_val

            if prot_id in modelable_prots:
                mutation_values[prot_id] += mut_val

        # keep track of all missense mutations in modelable proteins
        if mut_type in ['Missense'] and prot_id in modelable_prots:
            mutation_locations[prot_id].append((mut_pos, mut_val))

    return mutation_locations, mutation_values, total_mutational_value, total_mutations


####################################################################################################
# COVARIANCE COMPUTATIONS
####################################################################################################

def get_indices_from_interval(intervals_str):
    """
    :param intervals_str: comma-separated intervals (e.g., "2-10,20-45,90-100")
    :return: set of all indices as specified by the intervals string
    """

    current_indices = set()

    for current_interval in intervals_str.split(','):
        start_index, end_index = map(int, current_interval.split('-')[:2])
        for current_index in xrange(start_index, end_index + 1):
            current_indices.add(current_index)

    return current_indices


####################################################################################################

def covariance_to_correlation(cov_matrix):
    """
    :param cov_matrix: numpy array of shape NxN (matrix) of covariance values
    :return: numpy array of shape NxN (matrix) with the corresponding correlations
    """

    dim = len(cov_matrix)  # dimension of the square matrix
    std_deviations = [sqrt(cov_matrix[i][i]) for i in xrange(dim)]  # standard deviations from the diagonal

    # new blank correlation matrix to fill in:
    corr_matrix = np.zeros((dim, dim))  # default entry is np.float64(0)

    for i in xrange(dim):
        corr_matrix[i][i] = 1.
        for j in xrange(i + 1, dim):
            normalizer = std_deviations[i] * std_deviations[j]
            corr_matrix[i][j] = cov_matrix[i][j] / (normalizer if normalizer > 0 else 1)
            corr_matrix[j][i] = corr_matrix[i][j]

    return corr_matrix


####################################################################################################

def get_observed_mu_covariance(mutation_indices, weightfile, restriction='none'):
    """
    :param mutation_indices: list of mutated indices in this protein
    :param weightfile: full path to a file containing all binding potential weight tracks for the given protein
    :param restriction: str indicating a particular subset of tracks to consider (e.g., interaction, conservation, dom)
    :return: the expected score for each track, the observed scores for each track, the complete covariance
             matrix of all tracks with 1+ mutations, and the ordered list of track IDs that we have scores for
    """

    # Store information for all protein tracks with 1+ mutations
    covariances = {}  # track_id -> track_id -> score
    expected = {}  # track_id -> expected binding score
    observed = {}  # track_id -> observed binding score
    positive_mutations = {}  # track_id -> total mutations falling into positively-weighted positions
    intervals = {}  # track_id -> string of intervals (e.g., 2-10,20-45,90-100); needed if expected_muts == 'observed'
    track_names = {}  # track_id -> track_name

    # get list of domains that occur 40+ times in this protein (we will ONLY consider aggregate tracks in these cases)
    aggregate_names = None
    if restriction in ['none', 'interaction', 'nointeraction', 'domain', 'nodomain', 'noconservation', 'nowholegene',
                       'intercons', 'interdom', 'interwholegene', 'domcons', 'domwholegene']:
        aggregate_names = aggregate_domain_tracks(weightfile)

    weightfile_handle = gzip.open(weightfile) if weightfile.endswith('gz') else open(weightfile)
    for wline in weightfile_handle:
        if wline.startswith('#'):
            continue

        (track_id, track_name, intervals_yi, binding_weights,
         expected_yi, variance_yi, covariances_yizi) = wline[:-1].split('\t')[1:8]
        minimum_mut_count = max(int(wline[:-1].split('\t')[8]), 1) if len(wline[:-1].split('\t')) > 8 else 1

        # restrict to a subset of tracks based on the restriction passed in
        if not check_restrictions(track_name, restriction, aggregate_names):
            continue

        # skip tracks without the minimum required mutations for well-behaved Z-scores
        current_indices = get_indices_from_interval(intervals_yi)
        track_mutations = [(mut_index, mut_burden) for mut_index, mut_burden in mutation_indices
                           if mut_index in current_indices]
        del current_indices
        if len(track_mutations) < 1 or (len(track_mutations) < minimum_mut_count and
                                        track_name_to_classification(track_name) != 'conservation'):
            continue

        # observed binding score:
        weight_vector = {int(locweight.split(':')[0]): float(locweight.split(':')[1])
                         for locweight in binding_weights.split(',')}
        observed_val = sum([mut_burden * weight_vector.get(mut_index, 0.) for mut_index, mut_burden in track_mutations])

        if not observed_val > 0:
            continue

        # expected binding score:
        total_mutations = sum([mut_burden for _, mut_burden in track_mutations])
        expected_val = total_mutations * float(expected_yi)

        # note that we restrict to positive Z-scores only (for combination later)
        if not observed_val > expected_val:
            continue

        observed[track_id] = observed_val
        expected[track_id] = expected_val
        intervals[track_id] = intervals_yi  # interval that this track spans (for overlap calculations)
        track_names[track_id] = track_name  # name of track (to return name and classify track)
        # total value and number of positively weighted mutations (for confidence weights of tracks)
        positive_mutations[track_id] = (sum([mut_burden for mut_index, mut_burden in track_mutations
                                             if weight_vector.get(mut_index, 0.) > 0.]),
                                        len([mut_index for mut_index, _ in track_mutations
                                             if weight_vector.get(mut_index, 0.) > 0.]))

        # variance of binding score
        if track_id not in covariances:
            covariances[track_id] = {}
        squared_mutations = sum(
            [mut_burden * mut_burden for _, mut_burden in track_mutations])  # muts are not all = "1"
        covariances[track_id][track_id] = squared_mutations * float(variance_yi)

        # covariances between binding scores:
        if len(covariances_yizi) < 1:  # last track, no need to store its covariances to other tracks
            continue

        for other_id, cov in map(lambda x: x.split(':')[:2], covariances_yizi.split(',')):
            cov_score, likelihood_track, likelihood_other = map(float, cov.split('|')[:3])

            if other_id not in covariances:
                covariances[other_id] = {}

            covariances[track_id][other_id] = (cov_score, likelihood_track, likelihood_other)
            covariances[other_id][track_id] = (cov_score, likelihood_other, likelihood_track)

    weightfile_handle.close()

    # nothing passed! return empty:
    if len(observed.keys()) < 1:
        return {}, {}, np.empty(shape=(0, 0), dtype=float), [], {}, {}

    # finalize the covariance matrix
    sorted_tracks = [(observed[track_id] - expected[track_id], track_id) for track_id in observed.keys()]
    sorted_tracks.sort(reverse=True)  # sort by largest -> smallest effect size
    ids_to_check = [track_id for _, track_id in sorted_tracks]

    # add the very first track (variance and ID) to our matrix:
    covariance_matrix = [[covariances[ids_to_check[0]][ids_to_check[0]]]]
    final_ids = ids_to_check[:1]

    # only one track passed! no covariance matrix necessary:
    if len(ids_to_check) < 2:
        return expected, observed, np.array(covariance_matrix), final_ids, positive_mutations, track_names

    # incrementally add in other tracks if the resulting matrix remains positive definite
    for other_id in ids_to_check[1:]:
        update_covariance_matrix(covariance_matrix,  # current covariance matrix
                                 final_ids,  # list of tracks already included in covariance matrix
                                 other_id,  # new track id we would like to add
                                 covariances[other_id],  # how this new track covaries with all other tracks
                                 intervals,  # intervals that each track spans
                                 mutation_indices)  # indices where we had 1+ mutation

    return expected, observed, np.array(covariance_matrix), final_ids, positive_mutations, track_names


####################################################################################################

def update_covariance_matrix(current_cov, existing_tracks, other_id, covariances, intervals, mutation_indices):
    """
    NOTE:
    If the values are "too close" (i.e., covariance is very high?) between two tracks, then for all practical
    purposes (i.e., floating point precision), points from the two tracks are collinear. This means that the
    covariance matrix will not be positive DEFINITE, and thus we cannot take the Cholesky decomposition to
    speed up the process of finding the matrix inverse. The only possible fix is to check the addition of each
    track and make sure that the resulting matrix has only positive eigenvalues

    :param current_cov: list of N arrays of length N each corresponding to the existing (working) covariance matrix
    :param existing_tracks: ordered list of track ids corresponding to the tracks included in the existing cov. matrix
    :param other_id: track id of the new track that we are going to attempt to add
    :param covariances: dictionary of other_id -> (covariance score, likelihood of other_id, likelihood of current_id)
    :param intervals: dictionary of track_id -> string e.g., "1-10,20-25,90-120"
    :param mutation_indices: just incase "expected_muts" == 'observed', we'll need to know where mutations are falling
                             to determine the observed number of mutations
    :return: an updated covariance matrix and updated list of existing tracks, and a boolean for whether the addition
             was successful of not.
    """

    current_indices = get_indices_from_interval(intervals[other_id])  # indices for this new track
    current_muts = [(mut_index, mut_burden) for mut_index, mut_burden in mutation_indices
                    if mut_index in current_indices]  # subset of mutations falling into this new track

    new_track_row = []  # keep track of the newest (last) row to be added to the covariance matrix

    # for each existing track in the covariance matrix, determine if it has any overlap with this new track,
    # and include the covariance if so
    for i, track_id in enumerate(existing_tracks):

        # note that "covariances" contains how this current track ("other_id") covaries with all other tracks
        cov_score, likelihood_other, likelihood_track = covariances[track_id]

        # if the covariance or likelihoods are not all >0, then the overall covariance will of course be 0
        if not cov_score > 0 and likelihood_track > 0 and likelihood_other > 0:
            new_track_value = 0.

        else:  # otherwise, we compute covariance based on the number of mutations falling into the overlap region:
            overlap_indices = current_indices.intersection(get_indices_from_interval(intervals[track_id]))
            overlap_mutations = sum([mut_burden * mut_burden for mut_index, mut_burden in current_muts
                                     if mut_index in overlap_indices])

            new_track_value = cov_score * float(overlap_mutations)

        # extend the matrix in both directions: new row AND existing columns:
        new_track_row.append(new_track_value)
        current_cov[i].append(new_track_value)

    # finally, add the variance for this new track as the last element of the most recent row:
    new_track_row.append(covariances[other_id])  # already scaled by the total number of mutations in the track
    current_cov.append(new_track_row)  # matrix is now square again
    existing_tracks.append(other_id)  # note that we've added this newest track


####################################################################################################
# CALCULATE COMBINED Z-SCORE
####################################################################################################

def calculate_wholegene_zscore(gene_probability, current_mut_val, total_mut_cnt, total_mut_val):
    """
    :param gene_probability: relative likelihood of this gene harboring a nonsynonymous mutation
    :param current_mut_val: total sum of nonsynonymous mutation values landing in this gene
    :param total_mut_cnt: total sum of nonsynonymous mutation COUNTS landing across all modelable genes
    :param total_mut_val: total sum of nonsynonymous mutation VALUES landing across all modelable genes
    :return: z-score indicating the enrichment or depletion of mutations falling into this gene
    """

    # determine how to scale the total mutation count to get reasonably-sized Z-scores

    exp_mutval = total_mut_val / total_mut_cnt  # average mutation value
    scaled_mutcount = sqrt(total_mut_val) / exp_mutval  # scale the total mutations down to its square root
    total_mutval = scaled_mutcount * exp_mutval  # reset the total mutational value
    totsq_mutval = scaled_mutcount * exp_mutval * exp_mutval  # reset the total mutational squared values
    scale_factor = sqrt(total_mut_val) / total_mut_val  # scale down the eventual observed mut counts

    # whole gene track looks like: [0,0,0,0,...,1,...,0,0,0,0]. The sum of this AND the sum of these entries squared
    # are both 1, which is why the variance calculation works out below (the likelihood of landing on the 1 is same)

    mut_expected = gene_probability
    mut_variance = mut_expected - mut_expected ** 2

    prot_observed = current_mut_val * scale_factor
    prot_expected = total_mutval * mut_expected
    prot_variance = totsq_mutval * mut_variance

    return (prot_observed - prot_expected) / sqrt(prot_variance if prot_variance > 0 else 1.)


####################################################################################################

def protein_ztransform(mutation_indices, weightfile, current_mutational_value, total_mutational_value,
                       total_mutation_count, restriction='none'):
    """
    :param mutation_indices: list of mutated indices and their corresponding values in this protein
    :param weightfile: full path to a file containing all binding potential weight tracks for the given protein
    :param current_mutational_value: sum of mutation values falling into this particular protein
    :param total_mutational_value: total sum of mutation values falling across all modelable proteins
    :param total_mutation_count: total number of distinct mutational events occurring across all modelable proteins
    :param restriction: str indicating a particular subset of tracks to consider (e.g., interaction, conservation, dom)
    :return: a combined Z-score
    """

    final_zscores = []  # record the track names and Z-scores for eventual mechanism interpretation

    # (1) try to read from weightfile:
    if not weightfile:
        return 0., ';'.join(final_zscores)

    try:
        weightfile_handle = gzip.open(weightfile) if weightfile.endswith('gz') else open(weightfile)
        for _ in weightfile_handle:
            continue
        weightfile_handle.close()
    except (IOError, StopIteration) as _:  # weightfile is corrupted
        return 0., ';'.join(final_zscores)

    # (2) start with the whole gene Z-score if specified
    wholegene_zscore = 0.
    if restriction not in ['interaction', 'domain', 'conservation', 'interdom', 'intercons', 'domcons', 'nowholegene']:
        gene_probability = None
        weightfile_handle = gzip.open(weightfile) if weightfile.endswith('gz') else open(weightfile)
        for wt_line in weightfile_handle:
            if wt_line.startswith('# Relative Mutability & Total Genes Evaluated ='):
                gene_probability = float(wt_line.strip().split()[-2])
                break
        weightfile_handle.close()

        if gene_probability:
            wholegene_zscore = calculate_wholegene_zscore(gene_probability,
                                                          current_mutational_value,
                                                          total_mutation_count,
                                                          total_mutational_value)
            if wholegene_zscore > 0.:
                final_zscores = ['WholeGene_NatVar|' + str(wholegene_zscore)]

    # (3) now process all other tracks:
    (expected, observed, covariance_matrix,
     final_ids, positive_mutation_count, track_names) = get_observed_mu_covariance(mutation_indices,
                                                                                   weightfile,
                                                                                   restriction)
    # no additional tracks to include/integrate:
    if len(final_ids) < 1:
        if len(final_zscores) > 0:
            return wholegene_zscore, ';'.join(final_zscores)
        return 0., ';'.join(final_zscores)

    # give the interaction, domain, conservation, and natvar tracks relative weights of 1/4 each
    total_weight = {'interaction': 0., 'domain': 0.}
    track_classes = {}
    for track_id in final_ids:  # track ID -> (type, name)
        track_classes[track_id] = track_name_to_classification(track_names[track_id])
        if track_classes[track_id] in total_weight:
            current_track_value = positive_mutation_count[track_id][1]  # this should be non-zero by definition
            total_weight[track_classes[track_id]] += sqrt(current_track_value) if current_track_value > 0 else 0.

    # convert the covariance matrix to a correlation matrix:
    correlation_matrix = covariance_to_correlation(covariance_matrix)

    # get the vector of Z-scores and the vector of weights:
    # MC Whitlock (2005) "Combining probability from independent tests: the weighted Z-method is superior to
    # Fisher's approach" (http://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2005.00917.x/full)
    zscores, weights = [], []
    for i, track_id in enumerate(final_ids):

        std_dev = sqrt(covariance_matrix[i][i]) if covariance_matrix[i][i] > 0 else 1.
        unweighted_zscore = abs(observed[track_id] - expected[track_id]) / std_dev
        zscores.append(unweighted_zscore)

        # weight proportionally such that each class gets equal (i.e., sum to 1) weighting
        if track_classes[track_id] in total_weight:
            weights.append((sqrt(positive_mutation_count[track_id][1])
                            if positive_mutation_count[track_id][1] > 0 else 0.) /
                           total_weight[track_classes[track_id]])
        else:
            weights.append(1.)
        final_zscores.append(edit_track_name(track_names[track_id]) + '|' + str(unweighted_zscore))

    # combine the whole gene zscore, too:
    if restriction not in ['interaction', 'domain', 'conservation', 'interdom', 'intercons', 'domcons',
                           'nowholegene'] and wholegene_zscore and wholegene_zscore > 0.:
        # there is no within-gene variance for the whole gene track, so its covariance with all other tracks must be 0:
        correlation_matrix = np.append(np.vstack([correlation_matrix, [0.] * len(zscores)]),
                                       np.array([0.] * len(zscores) + [1.]).reshape((len(zscores) + 1, 1)), 1)
        zscores.append(wholegene_zscore)
        weights.append(1.)

    if not sum(weights) > 0:  # this only happens if nothing passed our filter (i.e., impossible)
        if len(final_zscores) > 0 and final_zscores[0].startswith('WholeGene_NatVar|'):
            return wholegene_zscore, ';'.join(final_zscores)
        return 0., ';'.join(final_zscores)  # nothing passed our weighting filter..

    # finally, compute the combined z-score for this protein:
    weights = np.array(weights)
    zscores = np.array(zscores)

    numerator = weights.dot(zscores)
    denominator = weights.dot(correlation_matrix).dot(weights.T)
    denominator = sqrt(denominator) if denominator > 0 else 1.
    combined_zscore = float(numerator / denominator)

    return combined_zscore, ';'.join(final_zscores)


####################################################################################################
# CREATE FINAL OUTPUT FILE
####################################################################################################

def mapping_gene_to_name(annotation_file):
    """
    :param annotation_file: full path to a tab-delimited file with Ensembl gene ID and primary gene names
    :return: dictionary Ensembl gene ID -> primary gene name
    """

    gene_to_name = {}
    annot_handle = gzip.open(annotation_file, 'rt') if annotation_file.endswith('gz') else open(annotation_file)
    header = None
    for annot_line in annot_handle:
        if annot_line.startswith('#'):
            continue
        elif not header:
            header = annot_line[:-1].split('\t')
            continue
        gene_zscore = annot_line[:-1].split('\t')
        gene_id = gene_zscore[header.index('ensembl_gene_id')]
        gene_name = gene_zscore[header.index('primary_gene_name(s)')]
        gene_to_name[gene_id] = gene_name
    annot_handle.close()

    return gene_to_name


####################################################################################################

def mapping_gene_to_driver(annotate_drivers, driver_annotation_file):
    """
    :param annotate_drivers: boolean indicating whether drivers should be annotated at all
    :param driver_annotation_file: full path to a tab-delimited file with Ensembl gene ID and driver status
    :return: dictionary Ensembl gene ID -> driver status AND list of driver status definitions (commented)
    """

    gene_to_cancer = {}
    cancer_header = []
    if annotate_drivers:
        if driver_annotation_file.endswith('gz'):
            annot_handle = gzip.open(driver_annotation_file, 'rt')
        else:
            annot_handle = open(driver_annotation_file)
        header = None
        for annot_line in annot_handle:
            if annot_line.startswith('##'):
                cancer_header.append(annot_line)
                continue
            elif annot_line.startswith('#'):
                continue
            elif not header:
                header = annot_line[:-1].split('\t')
                continue
            gene_zscore = annot_line[:-1].split('\t')
            gene_id = gene_zscore[header.index('ensembl_gene_id')]
            cancer_status = gene_zscore[header.index('cancer_driver_status')]
            gene_to_cancer[gene_id] = cancer_status
        annot_handle.close()

    return gene_to_cancer, cancer_header


####################################################################################################

def reformat_results(input_files, concatenated_output_file, annotation_file, annotate_drivers=False,
                     driver_annotation_file=None):
    """
    :param input_files: set of all files containing genes, their cancer status, and their scores
    :param concatenated_output_file: single output file containing the combined output from the input files
    :param annotation_file: full path to tab delimited list of gene ID, primary gene names
    :param annotate_drivers: boolean indicating whether cancer driver status should be included in output
    :param driver_annotation_file: full path to a tab-delimited file with Ensembl gene ID, driver status
    :return: None
    """

    # Get Ensembl gene ID -> (best isoform score, total time to run across all isoforms)
    gene_to_score = {}
    header = ['prot_id', 'gene_id', 'score', 'total_time', 'track_zscores']
    for infile in input_files:
        if not os.path.isfile(infile):
            continue

        infile_handle = gzip.open(infile) if infile.endswith('gz') else open(infile)
        for prot_line in infile_handle:
            if '\t' not in prot_line:
                continue
            v = prot_line[:-1].split('\t')

            # (1) get the gene name and the combined Z-score for this isoform
            gene_id = v[header.index('gene_id')]
            gene_zscore = float(v[header.index('score')])

            # (2) keep track of the time it took to run this particular protein:
            total_seconds = 0.
            if len(v) > header.index('total_time'):
                total_seconds = float(v[header.index('total_time')])

            # (3) and also which tracks had positive Z-scores:
            track_names = {}
            if len(v) > header.index('track_zscores') and len(v[header.index('track_zscores')].strip()) > 0:
                bd_tracks = v[header.index('track_zscores')].strip()
                track_names = {tname.split('|')[0]: float(tname.split('|')[1]) for tname in bd_tracks.split(';')}

            # (4) include this gene if it hasn't yet been observed:
            if gene_id not in gene_to_score:
                gene_to_score[gene_id] = {'score': gene_zscore,
                                          'time': total_seconds,
                                          'tracks': track_names}

            # (5) update the total seconds, score, and track names otherwise
            else:
                # (5a) update the score (i.e., take the maximum)
                gene_to_score[gene_id]['score'] = max(gene_to_score[gene_id]['score'], gene_zscore)

                # (5b) tack onto the total time
                gene_to_score[gene_id]['time'] += total_seconds

                # (5c) tack on the track names (and keep track of the maximum for corresponding tracks)
                for bd_track, zscore in track_names.items():
                    if bd_track not in gene_to_score[gene_id]['tracks']:
                        gene_to_score[gene_id]['tracks'][bd_track] = zscore
                    else:
                        gene_to_score[gene_id]['tracks'][bd_track] = max(gene_to_score[gene_id]['tracks'][bd_track],
                                                                         zscore)
        infile_handle.close()

    # (6) get the total time to run ALL genes
    total_runtime = sum([gv['time'] for gv in gene_to_score.values()])

    # (7) get the final list of sorted genes, reformatting the track names and time to be strings
    final_sorted_genes = sorted([(gene_vals['score'],
                                  gene_id,
                                  ';'.join([str(tn) + '|' + str(zs) for tn, zs in gene_vals['tracks'].items()]),
                                  reformat_time(gene_vals['time']))
                                 for gene_id, gene_vals in gene_to_score.items()], reverse=True)

    # (8) finally write out results!
    gene_to_name = mapping_gene_to_name(annotation_file)  # Ensembl gene ID -> primary gene name
    (gene_to_cancer,  # # Ensembl gene ID -> cancer driver status if desired
     cancer_header) = mapping_gene_to_driver(annotate_drivers, driver_annotation_file)

    concat_outhandle = open(concatenated_output_file, 'w')
    concat_outhandle.write('\n'.join(['# PertInInt (v0) results',
                                      '# Time to run = ' + reformat_time(total_runtime)] + cancer_header +
                                     ['\t'.join(['cancer_status', 'gene_name', 'score', 'runtime',
                                                 'per_track_scores'])]) + '\n')

    for gene_zscore, gene_id, track_zs, gene_runtime in final_sorted_genes:
        cancer_status = gene_to_cancer.get(gene_id, '')
        gene_full_name = gene_id + (',' if gene_id in gene_to_name else '') + gene_to_name.get(gene_id, '')
        concat_outhandle.write(
            '\t'.join([cancer_status, gene_full_name, str(gene_zscore), gene_runtime, track_zs]) + '\n')
    concat_outhandle.close()

    # (9) clean up the temporary (unformatted) output files
    for infile in input_files:
        if infile != concatenated_output_file:
            call(['rm', infile])


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    # ------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Return genes that are enriched for mutations in functional sites.')
    parser.add_argument('--maf_file', type=str, help='Mutation file (.maf format)', default=None)
    parser.add_argument('--out_file', type=str, help='Results output file', default=None)
    parser.add_argument('--ensembl_annotation_file', type=str, default='GRCh38_ensembl_gene_list.tsv.gz',
                        help='Tab-delimited list of Ensembl gene identifiers and their primary gene names')
    parser.add_argument('--track_path', type=str, help='Full path to directory containing track weight information',
                        default='track_weights/')

    parser.add_argument('--expression_file', type=str, default='TCGA_GRCh38_expressed-genes_TPM.tsv.gz',
                        help='Full path to a tab-delimited file containing lists of genes that are expressed in ' +
                             'particular tumor sample IDs')
    parser.add_argument('--driver_annotation_file', type=str, default='GRCh38_driver_gene_list.tsv.gz',
                        help='Tab-delimited list of Ensembl gene identifiers and their cancer driver status')

    parser.add_argument('--no_expression', dest='limit_expression', action='store_false', default=True,
                        help='Don\'t restrict mutations to those found in genes expressed at >0.1 TPM?')
    parser.add_argument('--no_driver_id', dest='annotate_drivers', action='store_false', default=True,
                        help='Don\'t annotate output with cancer driver gene status?')

    parser.add_argument('--restriction', type=str, help='Restrict to certain lines of functionality evidence',
                        default='none',
                        choices={'none', 'interaction', 'nointeraction', 'domain', 'nodomain',
                                 'conservation', 'noconservation', 'wholegene', 'nowholegene',
                                 'intercons', 'interdom', 'interwholegene', 'domcons', 'domwholegene', 'conswholegene'})
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------
    # (i) check if required files are all present
    if not args.maf_file or not os.path.isfile(args.maf_file):
        sys.stderr.write('* Could not open maf mutation file: ' + str(args.maf_file) + '\n' +
                         '* Please obtain a sample .maf file by running: \n' +
                         '    > if [ ! -d mafs ]; then mkdir mafs; fi\n' +
                         '    > AGGREGATE_CANCER=TCGA.Aggregate.muse.aggregated.somatic.maf.gz\n' +
                         '    > wget http://compbio.cs.princeton.edu/pertinint/$AGGREGATE_CANCER ' +
                         '-O mafs/$AGGREGATE_CANCER\n' +
                         '    > gzip -d mafs/$AGGREGATE_CANCER\n' +
                         '* Usage: python '+sys.argv[0]+' --maf_file <path_to_maf_file> ' +
                         '--out_file <path_to_output_file>\n' +
                         'Exiting.\n')
        sys.exit(1)

    if not os.path.isdir(args.track_path):
        sys.stderr.write('* Could not find directory to precomputed tracks: ' + args.track_path+'\n' +
                         '* Please obtain this directory by running: \n' +
                         '    > wget http://compbio.cs.princeton.edu/pertinint/PertInInt-tracks_v0.tar.gz\n' +
                         '    > tar -xvzf PertInInt-tracks_v0.tar.gz\n' +
                         '* Usage: python '+sys.argv[0]+' --maf_file <path_to_maf_file> ' +
                         '--out_file <path_to_output_file> --track_path <path_to_precomputed_tracks>\n' +
                         'Exiting.\n')
        sys.exit(1)

    if not os.path.isfile(args.ensembl_annotation_file):
        sys.stderr.write('* Could not open gene ID mapping file: '+args.ensembl_annotation_file + '\n' +
                         '* Please obtain this file by running: \n' +
                         '    > wget https://github.com/Singh-Lab/PertInInt/raw/master/' +
                         'GRCh38_ensembl_gene_list.tsv.gz\n' +
                         '* Usage: python ' + sys.argv[0] + ' --maf_file <path_to_maf_file> ' +
                         '--out_file <path_to_output_file> --expression_file <path_to_gene_mapping_file>\n' +
                         'Exiting.\n')
        sys.exit(1)

    if args.limit_expression and not os.path.isfile(args.expression_file):
        sys.stderr.write('* Could not open expression file: '+args.expression_file + '\n' +
                         '* Please obtain this file by running: \n' +
                         '    > wget http://compbio.cs.princeton.edu/pertinint/TCGA_GRCh38_expressed-genes_TPM.tsv.gz\n' +
                         '* Usage: python ' + sys.argv[0] + ' --maf_file <path_to_maf_file> ' +
                         '--out_file <path_to_output_file> --expression_file <path_to_expression_file>\n' +
                         'Proceeding WITHOUT restricting mutations to those falling into expressed genes.\n')
        args.limit_expression = False

    if args.annotate_drivers and not os.path.isfile(args.driver_annotation_file):
        sys.stderr.write('* Could not open driver annotation file: ' + args.driver_annotation_file + '\n' +
                         '* Please obtain this file by running: \n' +
                         '    > wget https://github.com/Singh-Lab/PertInInt/raw/master/GRCh38_driver_gene_list.tsv.gz\n' +
                         '* Usage: python ' + sys.argv[0] + ' --maf_file <path_to_maf_file> ' +
                         '--out_file <path_to_output_file> --driver_annotation_file <path_to_annotation_file>\n' +
                         'Proceeding WITHOUT annotating cancer driver genes in output file.\n')
        args.annotate_drivers = False

    # ------------------------------------------------------------------------------------------------
    # (ii) make sure that output file can be written to (exit quickly otherwise)
    for subdir in ['/'.join(args.out_file.split('/')[:ind]) for ind in xrange(2, args.out_file.count('/')+1)]:
        if not os.path.isdir(subdir):
            if call(['mkdir', subdir]):  # any code returned other than "0"
                sys.stderr.write('Could not write to '+args.out_file+'. Exiting\n')
                sys.exit(1)

    # ------------------------------------------------------------------------------------------------
    # (1) get paths to genes we can model
    sys.stderr.write('(1) Getting paths to proteins (in '+args.track_path+' directory) that can be modeled...\n')
    start = time.time()
    prot_to_trackfile, prot_to_geneid = track_weights_list(args.track_path)
    sys.stderr.write('    > finished in '+reformat_time(time.time()-start)+'\n')

    # ------------------------------------------------------------------------------------------------
    # (2) read in mutations for those genes that can modeled
    # NOTE: for whole gene tracks, we measure all nonsynonymous mutations, whereas for within-protein tracks,
    #       we are only interested in missense mutations
    sys.stderr.write('(2) Reading in mutation data from '+args.maf_file+'...\n' +
                     ('' if not args.limit_expression else
                      '    > Limiting to mutations in expressed genes found in: ' + args.expression_file + '\n'))
    start = time.time()

    (mut_locs,  # prot_id -> [(0-index position of missense mutation, mutational value), ...]
     mut_values,  # prot_id -> total mutational value (allelic fraction) for nonsynonymous mutations
     total_mut_value,  # total overall mutational value (allelic fraction) for nonsynonymous mutations
     total_mut_count) = process_mutations_from_maf(  # total overall mutation events (mutation count)
        args.maf_file,  # full path to .maf file
        set(prot_to_geneid.values()),  # set of Ensembl gene IDs with 1+ modelable proteins
        set(prot_to_trackfile.keys()),  # set of Ensembl protein IDs that can be modeled
        args.ensembl_annotation_file if args.limit_expression else None,
        args.expression_file if args.limit_expression else None  # path to expression file
    )
    sys.stderr.write('    > finished in '+reformat_time(time.time()-start)+'\n')

    # ------------------------------------------------------------------------------------------------
    # (3) start processing mutated proteins:
    sys.stderr.write('(3) Processing per-protein Z-scores for all mutated, modelable proteins...\n')
    start = time.time()

    out_handle = open(args.out_file+'-tmp', 'w')
    for mutated_protein, current_mutations in mut_locs.items():

        protein_start = time.time()  # start the clock to measure performance for this particular protein

        score, track_zscores = protein_ztransform(current_mutations,
                                                  prot_to_trackfile.get(mutated_protein, None),
                                                  mut_values.get(mutated_protein, 0.),
                                                  total_mut_value,
                                                  total_mut_count,
                                                  args.restriction)

        protein_total_time = time.time() - protein_start  # end clock to calculate total elapsed time (in seconds)

        # Write out these results:
        out_handle.write('\t'.join(map(str, [mutated_protein,
                                             prot_to_geneid.get(mutated_protein, ''),
                                             score,
                                             protein_total_time,
                                             track_zscores])) + '\n')
    out_handle.close()
    sys.stderr.write('    > finished in '+reformat_time(time.time()-start)+'\n')

    # ------------------------------------------------------------------------------------------------
    # (4) reformat per-protein results:
    sys.stderr.write('(4) Reformatting temporary output into '+args.out_file+'...\n' +
                     ('' if not args.annotate_drivers else
                      '    > Annotating with cancer driver status from: ' + args.driver_annotation_file + '\n'))
    start = time.time()
    reformat_results([args.out_file + '-tmp'], args.out_file, args.ensembl_annotation_file,
                     args.annotate_drivers, args.driver_annotation_file if args.annotate_drivers else None)
    sys.stderr.write('    > finished in '+reformat_time(time.time()-start)+'\n')
