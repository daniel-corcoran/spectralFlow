# Don't touch anything in this cell! Just run it and scroll to the bottom.

from math import sqrt
from scipy.spatial.distance import cosine
import os
import json
import operator
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

progress_dic = {} # Dictionary that tracks the percentage completion of certain jobs

def align_tokens(ta, tb, token_key):
    # Align two tokens to 0.002
    # First, detokenize and set asize mz pairs with delta <= 0.002
    dta = []  # mz
    daa = []  # abundance
    dtb = []
    dab = []
    for a, b, c in zip(ta, tb, token_key):
        if a > 0:
            dta.append(c)
            daa.append(a)
        if b > 0:
            dtb.append(c)
            dab.append(b)
    # Now we have a list of m/z and abundance..
    for x in range(len(dta)):
        for y in range(len(dtb)):
            if abs(dta[x] - dtb[y]) <= 0.002:
                # print(dta[x], dtb[y])
                dta[x] = dtb[y]
    my_token = []
    mz_ab_key = {}
    for mz, a in zip(dta, daa):
        mz_ab_key[mz] = a
    for t in token_key:
        if t in mz_ab_key:
            my_token.append(mz_ab_key[t])
        else:
            my_token.append(0)
    ta = my_token
    my_token = []
    mz_ab_key = {}
    for mz, a in zip(dtb, dab):
        mz_ab_key[mz] = a
    for t in token_key:
        if t in mz_ab_key:
            my_token.append(mz_ab_key[t])
        else:
            my_token.append(0)
    tb = my_token
    return ta, tb


def get_similarity(ta, tb, key):
    ta, tb = align_tokens(ta, tb, key)
    # score = 1 - cosine([sqrt(x) * mz for x, mz in zip(ta, key)], [sqrt(x) * mz for x, mz in zip(tb, key)])
    score = 1 - cosine([sqrt(x) for x in ta], [sqrt(x) for x in tb])
    return round(score, 5)


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def tokenize(sample, key):
    # Returns a token vector from the single sample.
    perfect = {}  # Values that immediately fit in the token vector.
    remainder = {}  # Values that don't immediately fit.
    sample_token = [0] * len(key)
    for mz, ab in zip(sample['mass'], sample['abundance']):
        if mz in key:  # Find the position and set it to abundance
            sample_token[key.index(mz)] = ab
        else:
            close = closest(key, mz)
            if abs(mz - close) <= 0.002:
                sample_token[key.index(close)] = ab
                # Here is where things get tricky.
            # If there is a value in the list within 0.002 of our value, we choose the closest one.
    return sample_token


def pretty_single_sample(single_sample, ab_cut):
    a_norm = max(single_sample['abundance'])
    # Normalize  masses to 572.2, highest RNA mass.
    single_sample['pepmass'] = single_sample['pepmass']
    single_sample['mass'] = [m for m in single_sample['mass']]
    single_sample['abundance'] = [m / a_norm for m in single_sample['abundance']]
    m = single_sample['mass']
    a = single_sample['abundance']

    single_sample['mass'] = []
    single_sample['abundance'] = []

    for i, j, in zip(m, a):
        if j >= ab_cut:
            single_sample['mass'].append(i)
            single_sample['abundance'].append(j)
    return single_sample


def yield_vector_from_file(path, ab_cut=0):
    with open(path) as source:
        idd = 0
        single_sample = {'mass': [], 'abundance': [], 'lines': [], 'checked': True, 'id': idd}  # One RT
        all_sample = []
        for line in source:
            single_sample['lines'].append(line)
            if 'END IONS' in line:  # Parse the current sample and clear cache
                single_sample = pretty_single_sample(single_sample, ab_cut)
                all_sample.append(single_sample)
                idd += 1
                single_sample = {'mass': [], 'abundance': [], 'lines': [], 'checked': True, 'id': idd}
            elif line[0] == 'T':
                single_sample['title'] = line.split()[0].replace('TITLE=', '')

                single_sample['scan'] = line.split()[-1][:-1].replace('scan=', '')
            elif line[0] in ['B', 'T', 'C', '\n']:  # This information is not used, so skip
                pass
            elif line[0] == 'R':
                single_sample['rt'] = float(line.split('=')[1])
            elif line[0] == 'P':
                single_sample['pepmass'] = float(line.split()[0][8:])
                # Only 8 decimals to filter out if abundance is counted
            else:
                mass, abundance = list(map(float, line.split()))
                if abundance >= ab_cut:
                    # print(mass, abundance)
                    single_sample['abundance'].append(round(abundance, 3))
                    single_sample['mass'].append(round(mass, 3))
        for z in all_sample:
            yield z


rmd = {}


def rna_mass(rna_string):
    # This function returns the average pep mass of an RNA modification
    # @rna: the name (string) of the RNA modification to query
    if rna_string in rmd:
        return rmd[rna_string]
    arr = []
    file = open("positive/" + rna_string + ".mgf", 'r')
    for line in file:
        if 'PEPMASS=' in line:  # Stop appending to the array when the ion stream ends
            mass = float(re.sub('\n', '', re.sub('PEPMASS=', '', line)))
            arr.append(mass)
    peps = np.asarray(arr)
    average = np.average(np.asarray(arr))
    rmd[rna_string] = average
    return average


def spectral_match_tkn(unknown_token, token_dic, mass_dic, pep_mass, mass_delta, check_mass):
    # Returns a dictionary of non-zero cosine similarities.
    rtrn_dic = {}
    non_key_tokens = {}
    for t in token_dic:
        if t != 'key':
            non_key_tokens[t] = token_dic[t]
    key = token_dic['key']

    for t in non_key_tokens:
        score = get_similarity(unknown_token, non_key_tokens[t], key)
        if str(score) != 'nan' and score != 0:
            if check_mass:
                if abs(mass_dic[t] - pep_mass) < mass_delta:
                    rtrn_dic[t] = score
            else:
                rtrn_dic[t] = score
    return rtrn_dic


def spectral_match(file_path='',
                   abundance_cutoff=0.0,
                   cosine_threshold=0.0,
                   mass_delta=0.001,
                   check_mass=True,
                   match_path='',
                   unmatched_path='',
                   ignore_modifications=[]):
    print(
        "Performing spectral search.\nFile: {}\nCosine Similarity Threshold: {}\nAbundance Cutoff: {}\nMass delta: {}".format(
            file_path, cosine_threshold, abundance_cutoff, mass_delta))
    count_dic = {}
    generator_object = yield_vector_from_file(file_path, ab_cut=abundance_cutoff)

    with open('json/token_{}.json'.format(abundance_cutoff)) as f:
        token_dic = json.load(f)
    key = token_dic['key']
    print('\nRT\tTop match\tCosine similarity')
    min_val = 1
    min_mod = ''

    match_write_string = ''
    unmatch_write_string = ''

    for x in generator_object:
        stats = spectral_match_tkn(tokenize(x, key), token_dic, x['pepmass'], mass_delta, check_mass)

        if len(stats):
            best_match = max(stats.items(), key=operator.itemgetter(1))[0]
            best_score = stats[best_match]
            if best_score >= cosine_threshold:
                if best_score < min_val:
                    min_val = best_score
                    min_mod = best_match
                if best_match not in count_dic:
                    count_dic[best_match] = 0
                count_dic[best_match] += 1

                if check_mass == False:
                    print(round(x['rt'] / 60, 2), '\t', best_match, '\t', best_score, '\t', x['scan'], '\t',
                          x['pepmass'])
                for line in x['lines']:
                    match_write_string += line
            else:
                # add to unmatched string
                for line in x['lines']:
                    unmatch_write_string += line
        else:
            # add to unmatched string
            for line in x['lines']:
                unmatch_write_string += line

    print("\n\nNumber of modifications detected: {}".format(len(count_dic)))
    for x in count_dic:
        print(x, count_dic[x])
    print("Lowest score: {} {}".format(min_val, min_mod))

    with open(match_path, 'w') as f:
        f.write(match_write_string)
    with open(unmatched_path, 'w') as f:
        f.write(unmatch_write_string)


def compare_mass(a, b, delta):
    if abs(a - b) <= delta:
        return True
    else:
        return False


def unique_ions(sample_a='',
                   sample_b='',
                   output_a='',
                   output_b='',
                   mass_delta=0.004):
    s_a = yield_vector_from_file(path=sample_a)
    s_b = yield_vector_from_file(path=sample_b)

    # Compare outputs from A False and B False
    # One output is the difference between A and B (What is in A that is not in B? (Only by pep mass))
    # Another output is same thing, but vice-versa
    # Keep the error of 0.004

    sample_a_unique = []
    sample_b_unique = []

    # Iterate through the mass lists and find values that aren't present in the other list
    for a in s_a:
        unique = True  # Assume the value is unique unless we can prove it isn't
        s_b = yield_vector_from_file(path=sample_b)
        for b in s_b:
            if compare_mass(a['pepmass'], b['pepmass'], mass_delta):
                unique = False
        if unique:
            sample_a_unique.append(a)

    s_b = yield_vector_from_file(path=sample_b)

    for b in s_b:
        unique = True  # Assume the value is unique unless we can prove it isn't
        s_a = yield_vector_from_file(path=sample_a)
        for a in s_a:
            if compare_mass(a['pepmass'], b['pepmass'], mass_delta):
                unique = False
        if unique:
            sample_b_unique.append(b)

    out_a_str = ''
    out_b_str = ''
    print("Unique samples in A")
    print("RT(min)\t", )
    for a in sample_a_unique:
        for line in a['lines']:
            out_a_str += line
    for b in sample_b_unique:
        for line in b['lines']:
            out_b_str += line
    with open(output_a, 'w') as f:
        f.write(out_a_str)
    with open(output_b, 'w') as f:
        f.write(out_b_str)


def mass_difference(file_path='',
                    matched_file = None,
                    unmatched_file = None,
                    threshold = 0.0,
                    mass=None,
                    mass_delta=None):
    # Find any samples that have peak m/zs  equal to (pep mass - mass)
    samples = yield_vector_from_file(path=file_path, ab_cut = threshold)
    relevant_samples = []
    irrelevant_samples = []
    relevant_mass = []
    for sample in samples:
        target = sample['pepmass'] - mass
        # We only care about the sample if there is a fragment with
        # mass 'target' in the mass vector
        is_sample_relevant = False
        for m in sample['mass']:
            if compare_mass(m, target, mass_delta):
                is_sample_relevant = True
                mass_a = m

        if is_sample_relevant:
            relevant_mass.append(mass)
            relevant_samples.append(sample)
        else:
            irrelevant_samples.append(sample)

    out_a_str = ''
    out_b_str = ''
    print("RT (min)\tPrecursor Mass\tIon Mass\tScan Number")
    for a, mass_a in zip(relevant_samples, relevant_mass):
        print(round(a['rt'] / 60, 2), a['pepmass'], mass, a['scan'])
        for line in a['lines']:
            out_a_str += line
    for b in irrelevant_samples:
        for line in b['lines']:
            out_b_str += line
    with open(matched_file, 'w') as f:
        f.write(out_a_str)
    with open(unmatched_file, 'w') as f:
        f.write(out_b_str)


def check_progress(ses_tkn):
    if ses_tkn in progress_dic:
        return progress_dic[ses_tkn]
    else:
        return 0


def web_api_hook(ses_tkn,
                 cos_threshold,
                 mass_check,
                 mass_delta,
                 abundance_cut,
                 include_reference,
                 include_sample):

    global progress_dic

    reference_gen = yield_vector_from_file(f'files/{ses_tkn}reference.mgf', ab_cut=abundance_cut)
    sample_gen = yield_vector_from_file(f'files/{ses_tkn}sample.mgf', ab_cut=abundance_cut)
    reference = [x for x in reference_gen if str(x['id']) in include_reference]
    sample = [x for x in sample_gen if str(x['id']) in include_sample]
    reference_id_title = {}
    for r in reference:
        reference_id_title[r['id']] = r['title']
    key = []

    for x in reference:
        for mz in x['mass']:
            if mz not in key:
                key.append(mz)

    # Now generate the token dic
    token_dic = {'key': key}
    count_dic = {}
    mass_dic = {}

    for x in reference:
        mass_dic[x['id']] = x['pepmass']
        token_dic[x['id']] = tokenize(x, key)

    print('\nRT\tTop match\tCosine similarity')
    min_val = 1
    min_mod = ''

    match_write_string = ''
    unmatch_write_string = ''

    results = []
    z = len(sample)
    for counter, x in enumerate(sample):
        progress_dic[ses_tkn] = counter / z
        stats = spectral_match_tkn(tokenize(x, key), token_dic, mass_dic, x['pepmass'], mass_delta, mass_check)

        if len(stats):
            #print(int(100 * check_progress(ses_tkn)), "%")
            best_match = max(stats.items(), key=operator.itemgetter(1))[0]
            best_score = stats[best_match]
            #print(stats)

            if best_score >= cos_threshold:

                results.append({'sample': x['scan'],
                                'reference': reference_id_title[best_match],
                                'cosine': round(best_score, 3),
                                'ab_cut': abundance_cut,
                                'sample_id': x['id'],
                                'reference_id': best_match,
                                'sample_rt_s': x['rt'],
                                'sample_rt_m': round(x['rt']/ 60, 2)})

                if best_score < min_val:
                    min_val = best_score
                    min_mod = best_match
                if best_match not in count_dic:
                    count_dic[best_match] = 0
                count_dic[best_match] += 1

    #print("\n\nNumber of modifications detected: {}".format(len(count_dic)))
    for x in count_dic:
        print(x, count_dic[x])
    #print("Lowest score: {} {}".format(min_val, min_mod))
    return results
