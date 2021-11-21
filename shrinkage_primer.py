import numpy as np
import random
import editdistance
import time

primers = []
NUM = 2000
DIM = 20
THRESH = 0.4 * DIM
CG_MIN = 0.45
CG_MAX = 0.55
TEMP_NUM = 20000
start_time = time.time()

while len(primers) < NUM:
    # randomly generate TEMP_NUM primers with DIM length and satisfy CG_MIN/MAX
    temp_primers = []
    while len(temp_primers) < TEMP_NUM:
        primer = ''.join(random.choice(
            ['A', 'C', 'G', 'T']) for _ in range(DIM))
        CG_COUNT = primer.count('C') + primer.count('G')
        if CG_COUNT < DIM * CG_MIN or CG_COUNT > DIM * CG_MAX:
            continue
        temp_primers.append(primer)

    # group primers with edit distance < 0.4L together
    group_primers = []
    while len(temp_primers) > 0:
        init_primer = temp_primers[0]
        group = []
        for primer in temp_primers:
            if editdistance.eval(init_primer, primer) < THRESH:
                group.append(primer)
        for primer in group:
            temp_primers.remove(primer)
        group_primers.append(group)

    # add primer to the final primers
    violate = False
    for group in group_primers:
        if len(primers) == 0:
            primers.append(group[0])
            continue
        for primer in primers:
            if editdistance.eval(group[0], primer) < THRESH:
                violate = True
                break
        if not violate:
            primers.append(group[0])
            if len(primers) >= NUM:
                break

    temp_time = time.time()
    print(len(primers), "primers have been found in",
          str(temp_time - start_time), "sec")

end_time = time.time()

# validation
has_wrong = False
for i in range(len(primers)):
    for j in range(len(primers)):
        if i == j:
            continue
        if editdistance.eval(primers[i], primers[j]) < THRESH:
            has_wrong = True
            break
    if has_wrong:
        break

if has_wrong:
    print("Validation fail!")
else:
    print("Validation succeed! Find all", len(primers),
          "primers in", str(end_time - start_time), "sec")
