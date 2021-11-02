from geneticalgorithm import geneticalgorithm as ga
import editdistance
import numpy as np

BEST = []
DIM = 20
THRESH = DIM*0.4
CG_MIN = 0.4 * DIM
CG_MAX = 0.55 * DIM


def feedback(candidate):
    cg_count = 0
    for base in candidate:
        if base == 1 or base == 2:
            cg_count += 1
    if cg_count < CG_MIN or CG_MAX < cg_count:
        return 0

    if len(BEST) == 0:
        BEST.append(candidate)
        return 0

    distances = [editdistance.eval(candidate, rep) for rep in BEST]
    all_exceed = all([dist > THRESH for dist in distances])

    if all_exceed:
        BEST.append(candidate)
        print(len(BEST))
        return 0
    total = sum(distances)
    cost = min(distances) - (1/total)
    return -cost


varbound = np.array([[1, 4]]*DIM)

model = ga(function=feedback, dimension=DIM, variable_type='int',
           variable_boundaries=varbound)

while len(BEST) < 2000:
    model.run()
    print(len(BEST))
