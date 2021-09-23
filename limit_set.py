import numpy as np
from numpy.linalg import inv
import random

# Return an array of points approximating the limit set of a group using a depth-first search.
#
# Arguments:
#  generators - list of 2x2 matrix generators
#  seed - complex points to map by the generators to produce the limit limit set
#  depth - maximal word length to generate
#  colored - if false, return the limit set as a list only; if true, return the limit set as a list of pairs (g,p) where p is a point and g is the initial letter of the word used to generate it
#            in the form of (the position in the generator array + 1) if g is a generator, and the negative of that if g is the inverse of a generator.
#  only_leaves - if true return images of seed under all the words, if false only under the longest words
def dfs(generators, seed, depth, coloured, only_leaves=False):
    # word_list is list of things of the form  (previous generator added, current word, initial generator decorator)
    decorated_gens = [(g, generators[g], g+1) for g in range(len(generators))]
    decorated_invs = [(g+len(generators), inv(generators[g]), -g-1) for g in range(len(generators))]
    word_list = [decorated_gens + decorated_invs]

    generators.extend([inv(g) for g in generators])

    for d in range(depth):
        print("Computing at depth " + str(d+1) + "/"+ str(depth))
        word_list.append([]);
        for old_word in word_list[d]:
            word_list[d+1].extend([(g, np.matmul(old_word[1], generators[g]), old_word[2]) for g in range(len(generators)) if g != old_word[0]])

    if only_leaves:
        word_list = word_list[-2:-1]

    seed = np.stack((seed,np.ones(len(seed))))
    word_list = flat_list = [item for sublist in word_list for item in sublist]
    limit_set_projective = [(np.matmul(w[1], seed), w[2]) for w in word_list]
    limit_set = [(p[0][0]/p[0][1], p[1]) for p in limit_set_projective]

    if coloured:
      return limit_set
    else:
      return limit_set[0]




# Return an array of points approximating the limit set of a group using a Markov chain search.
#
# Arguments:
#  generators - list of 2x2 matrix generators
#  seed - complex points to map by the generators to produce the limit limit set
#  depth - maximal word length to generate
#  colored - if false, return the limit set as a list only; if true, return the limit set as a list of pairs (g,p) where p is a point and g is the initial letter of the word used to generate it
#            in the form of (the position in the generator array + 1) if g is a generator, and the negative of that if g is the inverse of a generator.
#  reps - how many words to generate
def markov(generators, seed, depth, coloured, reps):
    decorated_gens = { g + 1: generators[g] for g in range(len(generators))}
    decorated_gens.update({-g - 1: inv(generators[g]) for g in range(len(generators))})
#    decorated_gens = { g + 1: generators[g] for g in range(len(generators))} | {-g - 1: inv(generators[g]) for g in range(len(generators))}

    seed = np.stack((seed,np.ones(len(seed))))
    limit_set_projective = [];

    for n in range(reps):
      first = 0
      word = np.identity(2)
      last = 0
      for d in range(depth):
          #print("Computing at depth " + str(d+1) + "/"+ str(depth) + " (iterate " + str(n+1) + "/" + str(reps) +")")
          admissable_keys = list(decorated_gens)
          try:
            admissable_keys.remove(-last)
          except:
            pass
          key = random.choice(admissable_keys)
          last = key
          if first == 0:
            first = key
          word = np.matmul(word,decorated_gens[key])
          limit_set_projective.append( (np.matmul(word, seed), first) )

    limit_set = [p[0]/p[1] for p in limit_set_projective]

    limit_set = [(p[0][0]/p[0][1], p[1]) for p in limit_set_projective]

    if coloured:
      return limit_set
    else:
      return limit_set[0]
