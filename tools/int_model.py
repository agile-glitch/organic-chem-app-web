"""The learned interference weights and the context-aware reactivity prior, as the ranker (tools/rel_fuse.py) and the page's exporter use them.

    IntModel('v1')    organic_chem_intmodel_v1.pkl  (tools/int_train.py)   score(flags, family) = the learned log-score of a site
    CtxModel('new')   organic_chem_ctx_new.pkl      (tools/int_context.py) features(...) = what a candidate consumed / left, judged with the competitors present
"""
import os
import pickle
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D

TMP = tempfile.gettempdir()
NSIG = 3
TRIM = 0.02


class IntModel:
    def __init__(self, tag='v1'):
        m = pickle.load(open(os.path.join(TMP, f'organic_chem_intmodel_{tag}.pkl'), 'rb'))
        self.keys, self.fams, w = m['keys'], m['families'], m['w']
        self.m = m
        nf = len(self.keys)
        self.nf = nf
        self.kidx = {k: i for i, k in enumerate(self.keys)}
        g = w[:nf + NSIG]
        self.g = g
        self.eff = {}
        self.d = {}
        for f, name in enumerate(self.fams):
            a = nf + NSIG + f * (nf + NSIG)
            self.d[name] = w[a:a + nf + NSIG]
            self.eff[name] = g + self.d[name]
        self.sigma = np.zeros((nf, NSIG))
        for s in D.SITE:
            if s['sigma'] is not None and s['key'] in self.kidx:
                self.sigma[self.kidx[s['key']], {'ortho': 0, 'meta': 1, 'para': 2}[s['key'].split('-')[1]]] = s['sigma']

    def total(self, key, family=None):
        """the learned effect of one pattern: its weight (global, or global + the family's) plus, for a ring substituent, rho(position) x sigma (the Hammett term)"""
        i = self.kidx[key]
        w = self.g if family is None else self.eff.get(family, self.g)
        t = float(w[i])
        sg = self.sigma[i]
        if sg.any(): t += float((w[self.nf:self.nf + NSIG] * sg).sum())
        return t

    def score_idx(self, idx, family):
        """idx = indices into self.keys of the flags"""
        if not len(idx): return 0.0
        w = self.eff.get(family, self.g)
        idx = np.asarray(idx, np.int64)
        return float(w[idx].sum() + (w[self.nf:self.nf + NSIG] * self.sigma[idx].sum(0)).sum())

    def score(self, flag_keys, family):
        return self.score_idx([self.kidx[k] for k in flag_keys if k in self.kidx], family)


class CtxModel:
    """p(group F consumed) = sigmoid(logit(the reagent's rate) + b[F] + beta[F] . (groups present) + mult[F] * (F present twice or more))"""
    def __init__(self, tag='new'):
        m = pickle.load(open(os.path.join(TMP, f'organic_chem_ctx_{tag}.pkl'), 'rb'))
        self.m = m
        self.fragments = m['fragments']
        # the page gets the model rounded to 4 decimals with coefficients under 0.02 dropped (tools/int_export.py export_ctx); the ranker is fitted with the same numbers
        beta = np.round(m['beta'].astype(np.float64), 4)
        beta[np.abs(beta) < TRIM] = 0.0
        self.b, self.beta, self.mult = np.round(m['b'].astype(np.float64), 4), beta, np.round(m['mult2'].astype(np.float64), 4)

    def rates(self, base_rates, s_counts):
        """the context-adjusted consumption probability of every fragment (the ones not present do not matter)"""
        s = np.asarray(s_counts)
        present = (s >= 1).astype(np.float64)
        off = np.log(base_rates / (1 - base_rates))
        z = off + self.b + self.beta @ present + self.mult * (s >= 2)
        return np.clip(1 / (1 + np.exp(-z)), 1e-3, 1 - 1e-3)

    def features(self, p, s_counts, c_counts):
        """(sum log-odds of the consumed groups, sum log(1 - p) of the groups present and left, consumed, formed) — tools/rel_reactivity.py Prior.features, with p from rates()"""
        s = np.asarray(s_counts); c = np.asarray(c_counts)
        present = s >= 1
        cons = s > c
        formed = c > s
        a = float(np.log(p[cons] / (1 - p[cons])).sum()) if cons.any() else 0.0
        b = float(np.log(1 - p[present & ~cons]).sum()) if (present & ~cons).any() else 0.0
        ll = float(np.log(p[cons]).sum() if cons.any() else 0.0) + b                     # the log-likelihood of what this candidate consumed and left, groups taken one by one
        return a, b, int(cons.sum()), int(formed.sum()), ll
