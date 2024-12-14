import json


def load_cats(dataset, case):
    with open('config/prior_cats.json') as f:
        cats = json.load(f)
    if (dataset == 'pbmc10k'):
        for i in range(4):
            cats[dataset][str(i)] = cats[dataset]['all'].copy()
    cats = cats[dataset][case]
    return cats

def f_beta_score(prc, rcl, beta=0.1):
    if prc + rcl == 0:
        return 0
    return (1 + beta**2) * (prc * rcl) / ((prc * beta**2) + rcl)
