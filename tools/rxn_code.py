"""The reaction code, version 2, in Python: the same definition as js/rxn_code.js (read that file for the format).

     reactants>agents>products | conc=37%w/w | cond=reflux, ethanol, 2 equiv

The first word is a reaction SMILES; a '>' in it separates the fields except in '->' (a dative bond, as RDKit reads it). The old
version 1 code, starting materials~reagent~conditions~products, is still read and is no longer written.

    parse(text) -> dict(ok, version, reactants, agents, products, cond, conc, extra, error)
    build(reactants='', agents='', products='', cond='', conc='') -> str   (molecule fields are already-joined SMILES lines)
    to_v2(text) -> the same reaction as a version 2 code (None if it cannot be read)
    fields(text) -> (reactants, agents, cond, products) as dot-joined strings, from either version (raises ValueError if unreadable)
    code_column(head) -> the reaction-code column of a learned_reactions.xlsx sheet, whichever version wrote the workbook
"""
import re

LEARNED_CODE_COL = 'Reaction code (starting materials>reagent>product | cond=conditions)'     # learned_reactions.xlsx (tools/learn_reactions.py)

_SEP = re.compile(r'(?<!-)>')            # a field separator, not the '>' of a dative bond '->'


def _split(s):
    return [x.strip() for x in str(s or '').split('.') if x.strip()]


def _tail(tail):
    cond, conc, extra = '', None, {}
    for seg in str(tail or '').split('|'):
        seg = seg.strip()
        i = seg.find('=')
        if i <= 0: continue
        k, v = seg[:i].strip().lower(), seg[i + 1:].strip()
        if k == 'cond': cond = v
        elif k == 'conc': conc = v
        else: extra[k] = v
    return cond, conc, extra


def parse(text):
    t = str(text or '').strip()
    if not t: return {'ok': False, 'error': 'empty'}
    first = t.split()[0]
    if not _SEP.search(first) and '~' in t:                          # version 1
        parts = [x.strip() for x in t.split('~')]
        if len(parts) > 4: return {'ok': False, 'error': 'A reaction code is starting materials~reagent~conditions~products (at most three ~)'}
        parts += [''] * (4 - len(parts))
        return {'ok': True, 'version': 1, 'reactants': _split(parts[0]), 'agents': _split(parts[1]), 'products': _split(parts[3]),
                'cond': parts[2], 'conc': None, 'extra': {}}
    rxn = _SEP.split(first)
    cond, conc, extra = _tail(t[len(first):])
    if len(rxn) == 1:
        return {'ok': True, 'version': 2, 'reactants': _split(first), 'agents': [], 'products': [], 'cond': cond, 'conc': conc, 'extra': extra}
    if len(rxn) != 3: return {'ok': False, 'error': 'A reaction SMILES has exactly two ">" (reactants>agents>products)'}
    return {'ok': True, 'version': 2, 'reactants': _split(rxn[0]), 'agents': _split(rxn[1]), 'products': _split(rxn[2]),
            'cond': cond, 'conc': conc, 'extra': extra}


def build(reactants='', agents='', products='', cond='', conc=''):
    rxn = f'{reactants or ""}>{agents or ""}>{products or ""}'
    if rxn == '>>' and not cond and not conc: return ''
    tail = (['conc=' + conc] if conc else []) + (['cond=' + cond] if cond else [])
    return rxn + (' | ' + ' | '.join(tail) if tail else '')


def to_v2(text):
    p = parse(text)
    if not p['ok']: return None
    return build('.'.join(p['reactants']), '.'.join(p['agents']), '.'.join(p['products']), p['cond'], p['conc'])


def code_column(head):
    return next(h for h in head if isinstance(h, str) and h.startswith('Reaction code'))


def fields(text):
    p = parse(text)
    if not p['ok']: raise ValueError(f'not a reaction code: {text!r} ({p["error"]})')
    return '.'.join(p['reactants']), '.'.join(p['agents']), p['cond'], '.'.join(p['products'])
