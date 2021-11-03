""" Generic methods for Kleinian groups.
"""

from multiprocessing import Pool
from mpmath import mp, fp
import numpy as np
import random
import functools
import gc

def _fast_inv(mat):
    """ Invert a 2x2 matrix *assuming it has det 1*.
    """
    assert((mat.rows,mat.cols) == (2,2))
    return mp.matrix([[mat[1,1],-mat[0,1]],[-mat[1,0],mat[0,0]]])

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
            word_list[d+1].extend([(g, old_word[1] * generators[g], old_word[2]) for g in range(len(generators)) if g != old_word[0]])

    if only_leaves:
        word_list = word_list[-2:-1]

    seed.rows = 2
    for i in range(seed.cols):
      seed[1,i] = 1
    word_list = [item for sublist in word_list for item in sublist]
    limit_set_projective = [(w[1] * seed, w[2]) for w in word_list]

    if coloured:
        return [([q[0]/q[1] for q in p[0].transpose()], p[1]) for p in limit_set_projective]
    else:
        return [[q[0]/q[1] for q in p[0].transpose()] for p in limit_set_projective]



def _dynamics_of_one_word(decorated_gens, seed, depth,rep):
    """ Generate a word of length depth using a Markov chain and return the orbits of seed under that word.

        Arguments and output format optimised for use in limit_set_markov not in user code.

        Arguments:
          decorated_gens -- list of pairs (n,X) where X is the nth generator (if n > 0) or the inverse of the nth generator (if n < 0).
          seed - complex points to map by the generators to produce the limit limit set, as points of CP.
          depth -- length of word to generate.
          rep -- this is the rep'th repetition in the loop.

        Returns:
          list of pairs (point,gen) where point is a complex point in the *affine* limit set and gen is the first letter of the word we generated
    """

    #print(f'_dynamics_of_one_word {rep}',flush=True)

    # We have to convert backwards and forwards between mpmath and numpy in order to cross the multiprocessing boundary,
    # due to issues like https://github.com/uqfoundation/dill/issues/238 and https://github.com/sympy/sympy/issues/11999.
    decorated_gens = {g: mp.matrix(m) for g,m in decorated_gens.items()}
    seed = mp.matrix(seed)

    random.seed()
    key = random.choice(list(decorated_gens))
    previous_letter = key
    image = decorated_gens[key] * seed

    orbit = [(p[0]/p[1], key) for p in image.T.tolist()]

    for d in range(1,depth):
        admissable_keys = list(decorated_gens)
        admissable_keys.remove(-previous_letter)
        key = random.choice(admissable_keys)
        previous_letter = key
        image = decorated_gens[key] * image
        orbit.extend([(np.clongdouble(p[0]/p[1]), key) for p in image.T.tolist()])

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
    decorated_gens = { (g + 1): generators[g] for g in range(len(generators))}
    decorated_gens.update({(-g - 1): _fast_inv(generators[g]) for g in range(len(generators))})


    # We have to convert backwards and forwards between mpmath and numpy in order to cross the multiprocessing boundary,
    # due to issues like https://github.com/uqfoundation/dill/issues/238 and https://github.com/sympy/sympy/issues/11999.
    decorated_gens = {g: np.matrix(fp.matrix(m).tolist(),dtype=np.clongdouble) for g,m in decorated_gens.items()}
    seed = np.array(seed,dtype=np.clongdouble)

    seed = np.stack((seed,np.ones(len(seed))))

    limit_set = []
    with Pool() as pool:
        for orbit in pool.imap(functools.partial(_dynamics_of_one_word, decorated_gens, seed, depth), range(reps)):
            for pair in orbit:
                yield (mp.mpmathify(pair[0]),pair[1])
            del orbit
    #for orbit in map(functools.partial(_dynamics_of_one_word, decorated_gens, seed, depth), range(reps)):
        #for pair in orbit:
            #yield pair
