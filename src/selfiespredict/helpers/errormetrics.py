import itertools

def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk

def topN_accuracy(n, PATH, TRUE_PATH):
    """Helper function that returns top n accuracy
    """
    with open(PATH) as f, open(TRUE_PATH,'r') as g:
        lines = f.read().splitlines()
        lines_true = g.read().splitlines()
        if len(lines) != len(lines_true)*n:
            raise ValueError()
    indices = []
    for index,i in enumerate(zip(grouper(n,lines),lines_true)):
        if i[1] in i[0]:
            indices.append(index)
    return (len(indices)/len(lines_true), indices)
