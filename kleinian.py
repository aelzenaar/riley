""" Generic methods for Kleinian groups.
"""

import numpy as np

import random
from multiprocessing import Pool
import functools
import gc

def _fast_inv(mat):
    """ Invert a 2x2 matrix *assuming it has det 1*.
    """
    assert(mat.shape == (2,2))
    return [[mat[1][1],-mat[0][1]],[-mat[1][0],mat[0][0]]]

def limit_set_dfs(generators, seed, depth, coloured, only_leaves=False):
    """ Return an array of points approximating the limit set of a group using a depth-first search.

        Arguments:
          generators - list of 2x2 matrix generators
          seed - complex points to map by the generators to produce the limit limit set
          depth - maximal word length to generate
          coloured - if false, return the limit set as a list only; if true, return the limit set as a list of pairs (g,p) where p is a point and g is the initial letter of the word used to generate it
                     in the form of (the position in the generator array + 1) if g is a generator, and the negative of that if g is the inverse of a generator.
          only_leaves - if True then return images of seed under all the words, otherwise return the images only under the longest words (default False)
    """

    # word_list is list of things of the form  (previous generator added, current word, initial generator decorator)
    decorated_gens = [(g, generators[g], g+1) for g in range(len(generators))]
    decorated_invs = [(g+len(generators), _fast_inv(generators[g]), -g-1) for g in range(len(generators))]
    word_list = [decorated_gens + decorated_invs]

    generators.extend([_fast_inv(g) for g in generators])

    for d in range(depth):
        print("Computing at depth " + str(d+1) + "/"+ str(depth))
        word_list.append([])
        for old_word in word_list[d]:
            word_list[d+1].extend([(g, np.matmul(old_word[1], generators[g]), old_word[2]) for g in range(len(generators)) if g != old_word[0]])

    if only_leaves:
        word_list = word_list[-2:-1]

    seed = np.stack((seed,np.ones(len(seed))))
    word_list = [item for sublist in word_list for item in sublist]
    limit_set_projective = [(np.matmul(w[1], seed), w[2]) for w in word_list]

    if coloured:
        return [([q[0]/q[1] for q in p[0].transpose()], p[1]) for p in limit_set_projective]
    else:
        return [[q[0]/q[1] for q in p[0].transpose()] for p in limit_set_projective]



def _dynamics_of_one_word(decorated_gens, seed, depth,rep):
    """ Generate a word of length depth using a Markov chain and return the orbits of seed under that word.

        Arguments and output format optimised for use in limit_set_markov not in user code.

        Arguments:
          depth -- length of word to generate
          _ -- ignored (means we can curry this function inside multiprocessing.Pool.map())

        Returns:
          list of pairs (point,gen) where point is a complex point in the *affine* limit set and gen is the first letter of the word we generated
    """

    print(f'_dynamics_of_one_word {rep}',flush=True)

    random.seed()
    key = random.choice(list(decorated_gens))
    first_letter = key
    previous_letter = key
    image = np.matmul(decorated_gens[key], seed)

    orbit = [(p[0]/p[1], first_letter) for p in image.transpose()]

    for d in range(1,depth):
        admissable_keys = list(decorated_gens)
        admissable_keys.remove(-previous_letter)
        key = random.choice(admissable_keys)
        previous_letter = key
        image = np.matmul(decorated_gens[key],image)
        orbit.extend([(p[0]/p[1], first_letter) for p in image.transpose()])
        del admissable_keys

    return orbit

def limit_set_markov(generators, seed, depth, reps):
    """ An iterator yielding points points approximating the limit set of a group using a Markov chain search.

        The number of points generated will be depth*reps. Each point yelded is a list of pairs (g,p) where p is a complex-valued limit
        set point and g is the initial letter of the word used to generate it in the form of (the position in the generator array + 1)
        if g is a generator, and the negative of that if g is the inverse of a generator.

        Arguments:
          generators - list of 2x2 matrix generators *with det 1*.
          seed - complex points to map by the generators to produce the limit limit set
          depth - maximal word length to generate
          reps - how many words to generate
    """
    decorated_gens = { g + 1: generators[g] for g in range(len(generators))}
    decorated_gens.update({-g - 1: _fast_inv(generators[g]) for g in range(len(generators))})

    seed = np.stack((seed,np.ones(len(seed))))

    limit_set = []
    with Pool(maxtasksperchild=100) as pool:
        for orbit in pool.imap(functools.partial(_dynamics_of_one_word, decorated_gens, seed, depth), range(reps)):
            for pair in orbit:
                yield pair
    #for orbit in map(functools.partial(_dynamics_of_one_word, decorated_gens, seed, depth), range(reps)):
        #for pair in orbit:
            #yield pair
