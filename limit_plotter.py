import numpy as np
import kleinian
import matplotlib.pyplot as plt
import readline

print("Limit set for Kleinian group")
gen_count = int(input("Number of generators: "))
generators = []
for i in range(gen_count):
    print("Enter matrix entries for generator " + str(i+1) + "/" + str(gen_count) + ".")
    print("E.g. the identity matrix should be entered as `1,0,0,1'.")
    generators.append(np.array([complex(t) for t in input("Entries: ").replace(" ","").split(",")]).reshape(2,2))

print("Which limit set algorithm? Choices: `dfs', `markov'")
algo = input("Algorithm: ")

depth = int(input("Maximal word length: "))
coloured = (input("Coloured output (y/n):") == "y")


if algo == "dfs":
    only_leaves = (input("Show just leaf points (y/n): ") == "y")
    ls = kleinian.limit_set_dfs(generators,np.array([1]),depth,coloured,only_leaves)
elif algo == "markov":
    reps = int(input("Number of points to generate: "))
    ls = kleinian.limit_set_markov(generators,np.array([1]),depth,coloured,reps)
else:
    raise RuntimeError

outtypes = input("Output types (f - file; s - onscreen display): ")
if 'f' in outtypes:
    filename = input("File name: ")


ls = [t for t in ls if np.real(t[0]) > -4 and np.real(t[0]) < 4]


if coloured:
    colours = {-2: 'r', -1:'b', 1:'g', 2:'y'}
    plt.scatter([np.real(t[0]) for t in ls],[np.imag(t[0]) for t in ls],c=[colours[t[1]] for t in ls],marker=".",s=.1,linewidths=0)
else:
    plt.scatter([np.real(t[0]) for t in ls],[np.imag(t[0]) for t in ls],marker=".",s=.1,linewidths=0)

plt.axis('equal')
plt.tight_layout()

if 'f' in outtypes:
    plt.savefig(filename,dpi=500)

if 's' in outtypes:
    plt.show()
