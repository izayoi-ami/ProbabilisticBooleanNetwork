def ddict():
    from collections import defaultdict
    return defaultdict(RDF)

def mdict():
    from collections import defaultdict
    return defaultdict(ddict)

def hamming_distance(x,y):
    return x.__xor__(y).digits(2).count(1)

