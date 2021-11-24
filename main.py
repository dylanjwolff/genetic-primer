import time
from geneticalgorithm import geneticalgorithm as ga
import editdistance
import numpy as np
import sys

BEST = []
DIM = 10
THRESH = 5
CG_MIN = 0.45 * DIM
CG_MAX = 0.55 * DIM
SKIP_COUNT = 0
EXEC_COUNT = 0
PRIMERS = 1000


def feedback(candidate):
    global EXEC_COUNT
    global SKIP_COUNT
    global BEST
    EXEC_COUNT += 1

    cg_count = 0
    for base in candidate:
        if base == 1 or base == 2:
            cg_count += 1
    if cg_count < CG_MIN or CG_MAX < cg_count:
        SKIP_COUNT += 1
        return 0

    if len(BEST) == 0:
        BEST.append(candidate)
        return 0

    distances = [editdistance.eval(candidate, rep) for rep in BEST]
    all_exceed = all([dist >= THRESH for dist in distances])

    if all_exceed:
        BEST.append(candidate)
        print(len(BEST))
        return 0
    total = sum(distances)
    cost = min(distances) - (1/total)
    return -cost


def check_cg_ok(candidate):
    global SKIP_COUNT
    cg_count = 0
    for base in candidate:
        if base == 1 or base == 2:
            cg_count += 1
    if cg_count < CG_MIN or CG_MAX < cg_count:
        SKIP_COUNT += 1
        return False

    return True


def squish_feedback(candidate):
    separated = np.reshape(candidate, (int(len(candidate)/DIM), DIM))

    distance = 0
    for row_num, row in enumerate(separated):
        if not check_cg_ok(row):
            return np.finfo('d').max

        for other_row in separated[row_num+1:]:
            this_distance = editdistance.eval(row, other_row)
            if this_distance < THRESH:
                return np.finfo('d').max
            distance += this_distance
    print(f"FOUND one {distance}")

    return distance


varbound = np.array([[1, 4]]*DIM)

search_model = ga(function=feedback, dimension=DIM, variable_type='int',
                  variable_boundaries=varbound)

while len(BEST) < PRIMERS:
    # while False:
    start = time.time()
    search_model.run()
    end = time.time()
    print(f"{BEST[-1]} candidate")
    print(f"{len(BEST)} candidates found")
    print(f"{end - start} elapsed")
    print(f"{SKIP_COUNT} of {EXEC_COUNT} skipped")
    EXEC_COUNT = 0
    SKIP_COUNT = 0
    current_primers = np.reshape(BEST, DIM*len(BEST))
    unsquished = squish_feedback(current_primers)

    varbound = np.array([[1, 4]]*DIM*len(BEST))
    squish_model = ga(function=squish_feedback, dimension=DIM*len(BEST), variable_type='int',
                      variable_boundaries=varbound)
    print(f"Unsquished total dist {unsquished}")
    start = time.time()
    squish_model.run([current_primers])
    end = time.time()
    print(f"{end - start} squish time elapsed")

# ina = np.array([0.]*40)
# r = squish_feedback(ina)
# print(r)
