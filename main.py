from threading import Timer
import time
from geneticalgorithm import geneticalgorithm as ga
import editdistance
import numpy as np
import sys
import os
from math import comb

LOGFILE_NAME = "log.txt"
BEST = []
DIM = 20 
THRESH = 9 # DIM*0.4
CG_MIN = 0.45 * DIM
CG_MAX = 0.55 * DIM
SKIP_COUNT = 0
EXEC_COUNT = 0
PRIMERS = 2000
MAX_TANGENT = 1
ADDED_THIS_ROUND = False


def feedback(candidate):
    global EXEC_COUNT
    global MAX_TANGENT
    global ADDED_THIS_ROUND
    global SKIP_COUNT
    global BEST
    EXEC_COUNT += 1

    if not check_cg_ok(candidate):
        SKIP_COUNT += 1
        return 0

    if len(BEST) == 0:
        ADDED_THIS_ROUND = True
        BEST.append(candidate)
        return 0

    distances = [editdistance.eval(candidate, rep) for rep in BEST]
    all_exceed = all([dist >= THRESH for dist in distances])

    if all_exceed:
        num_tangent = sum([dist == THRESH for dist in distances])
        if num_tangent > MAX_TANGENT:
            MAX_TANGENT += 1
            ADDED_THIS_ROUND = True
            BEST.append(candidate)
            return 0
        if num_tangent == len(BEST):
            ADDED_THIS_ROUND = True
            BEST.append(candidate)
            return 0
        cost = THRESH + 1 + num_tangent
        return -cost

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

class RepeatedTimer(object):
    def __init__(self, interval, function, *args, **kwargs):
        self._timer = None
        self.interval = interval
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.function(*self.args, **self.kwargs)
        # Starting timer only after fn has completed; may result in longer intervals than
        # specified, but no risk of spawning too many processes if function is slow.
        self.start()

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False

def log():
    global LOGFILE_NAME
    timestamp = time.time()
    primers = BEST.copy()
    num_primers = len(primers)

    # This is very heavy-weight, can comment out if we want smaller logging intervals
    # However, computation is done out-of-band on another thread so it shouldn't
    # impact performace too much and I thought these numbers might be interesting to
    # look at.
    # Feel free to add any other data that might be interesting
    total_dist = 0
    total_num_tang = 0
    for i, pa in enumerate(primers):
        for pb in primers[i+1:]:
            this_dist = editdistance.eval(pa, pb)
            total_dist += this_dist
            if this_dist == THRESH:
                total_num_tang += 1

    num_pairs = comb(num_primers, 2)
    avg_distance = total_dist / num_pairs
    avg_num_tang = total_num_tang / num_pairs

    # print(f"Logging calc took {time.time() - timestamp}")
    with open(LOGFILE_NAME, "a") as f:
        f.write(f"{timestamp},{num_primers},{avg_distance},{avg_num_tang}\n")


os.remove(LOGFILE_NAME)

varbound = np.array([[1, 4]]*DIM)

search_model = ga(function=feedback, dimension=DIM, variable_type='int',
           variable_boundaries=varbound)

RepeatedTimer(1, log)
while len(BEST) < PRIMERS:
    start = time.time()
    search_model.run()
    end = time.time()

    if not ADDED_THIS_ROUND:
        best_effort_cost = search_model.output_dict["function"]
        best_effort_primer = search_model.output_dict["variable"]
        num_tang = -(best_effort_cost + THRESH + 1)
        if best_effort_cost < -(THRESH + 1):
            BEST.append(best_effort_primer)
            MAX_TANGENT = num_tang
            print(f"Added Best effort! num tang {num_tang}");


    print(f"{end - start} elapsed")
    print(f"Found {len(BEST)} primers")
    print(f"{SKIP_COUNT} of {EXEC_COUNT} skipped")
    EXEC_COUNT = 0
    SKIP_COUNT = 0
    ADDED_THIS_ROUND = False
