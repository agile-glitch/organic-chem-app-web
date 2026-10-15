"""Keeps docs/AGILES1_format.md in step with the AGILES code (js/mol_code.js) and its reader (tools/agiles_read.py).

  python tools/agiles_doc.py            check: exits 1 and says what is out of date
  python tools/agiles_doc.py --write [--example-file F]   rewrite the document's reader copy and worked example
                                        (from the example code already in the document, or the code in file F)

Checked: the reader copy in the document is tools/agiles_read.py byte for byte; the worked example's printed output is
what the reader prints for the example code; every item name, label mode and identifier js/mol_code.js writes, and its
code prefix and format line, are described in the document. Run after every change to the AGILES code (CLAUDE.md).
"""
import io, os, re, sys
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
DOC = os.path.join(APP, 'docs', 'AGILES1_format.md')
READER = os.path.join(HERE, 'agiles_read.py')
JS = os.path.join(APP, 'js', 'mol_code.js')


def block(text, name):
    m = re.search(rf'<!-- {name}:begin -->\n(.*?)\n<!-- {name}:end -->', text, re.S)
    return m


def reader_output(code):
    sys.path.insert(0, HERE)
    import agiles_read as A
    buf = io.StringIO()
    with redirect_stdout(buf): A.show(A.parse(A.unpack(code)))
    return buf.getvalue().rstrip('\n')


def example_code(doc):
    m = block(doc, 'example')
    c = re.search(r'```\n((?:AGILES1|ACM1):[A-Za-z0-9_-]+)\n```', m.group(1)) if m else None
    return c.group(1) if c else None


def render(doc, code):
    reader = open(READER, encoding='utf-8').read().rstrip('\n')
    ex = f'```\n{code}\n```\n\nWhat `python agiles_read.py` prints for it:\n\n```\n{reader_output(code)}\n```'
    doc = re.sub(r'(<!-- example:begin -->\n).*?(\n<!-- example:end -->)', lambda m: m.group(1) + ex + m.group(2), doc, flags=re.S)
    doc = re.sub(r'(<!-- reader:begin -->\n).*?(\n<!-- reader:end -->)', lambda m: m.group(1) + f'```python\n{reader}\n```' + m.group(2), doc, flags=re.S)
    return doc


def problems(doc):
    out = []
    code = example_code(doc)
    if not code: out.append('the worked example has no code')
    elif render(doc, code) != doc: out.append('the reader copy or the worked example is not what tools/agiles_read.py gives (run --write)')
    js = open(JS, encoding='utf-8').read()
    names = set(re.findall(r"item\(NS \+ '([a-z_.]+)'", js))
    ids = re.search(r"const IDS = \[([^\]]*)\]", js).group(1); modes = re.search(r"const MODES = \[([^\]]*)\]", js).group(1)
    names |= {x.strip(" '") for x in ids.split(',')} | {'labels.' + x.strip(" '") for x in modes.split(',')}
    names.discard('labels.')
    for n in sorted(names):
        if f'AGILES.{n}' not in doc and f'`{n.replace("labels.", "")}`' not in doc: out.append(f'item AGILES.{n} (js/mol_code.js) is not described')
    prefix = re.search(r"const PREFIX = '([^']+)'", js).group(1); fmt = re.search(r"FORMAT = '([^']+)'", js).group(1)
    if prefix not in doc: out.append(f'the code prefix {prefix} is not in the document')
    if f'`{fmt}`' not in doc: out.append(f'the format line "{fmt}" is not in the document')
    return out


def main():
    doc = open(DOC, encoding='utf-8').read()
    if '--write' in sys.argv:
        i = sys.argv.index('--example-file') if '--example-file' in sys.argv else -1
        code = open(sys.argv[i + 1]).read().strip() if i >= 0 else example_code(doc)
        doc = render(doc, code)
        open(DOC, 'w', encoding='utf-8').write(doc)
    p = problems(doc)
    for x in p: print('OUT OF DATE:', x)
    print('docs/AGILES1_format.md is in step with js/mol_code.js and tools/agiles_read.py' if not p else f'{len(p)} problem(s)')
    sys.exit(1 if p else 0)


if __name__ == '__main__':
    main()
