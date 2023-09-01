
from itertools import groupby

def get_barriers_names(barriers_mod: dict) -> dict:
    '''Each barrier is divided into sections with constant diameter. Most barriers do of course have only one section, but in the preprocessing
        the barrier diameter is calculated based on casing etc and not given explicitely. Also top and base is divided into these sections in the class object.
        Original format:   barr1_0: barr1, barr1_1: barr1, barr2_0: barr2    where barr1 here is divided into the two sections barr1_0 and barr1_1
        New format         barr1: [barr1_0, barr1_1], barr2: [barr2_0]'''
    names = barriers_mod["b_name"]
    bnames = {key: [j for j, _ in list(value)] for key, value in groupby(names.items(), lambda x: x[1])}

    return bnames
