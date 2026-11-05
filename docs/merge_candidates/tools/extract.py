#!/usr/bin/env python3
"""Extract function definitions from NERDSS C++ sources (rough brace matcher).

Writes $DUPSCAN_OUT/funcs.json (default tools/out/): list of
  {file, name, qual, start, end, nlines, body_lines (comment-stripped, ws-normalized),
   tokens (identifier/number-normalized)}
"""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import json, os, re, sys


DIRS = ["src", "EXEs", "include"]
OUT = _os.path.join(OUTDIR, "funcs.json")

KEYWORDS = set("""if for while switch return sizeof catch else do new delete throw case
default static_cast dynamic_cast reinterpret_cast const_cast decltype alignof typeid""".split())


def strip_comments(src):
    """Replace comments with spaces (keep newlines), and string/char literals with "S"/'C' (same length newlines kept)."""
    out = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c == '/' and i + 1 < n and src[i + 1] == '/':
            j = src.find('\n', i)
            if j < 0: j = n
            out.append(' ' * (j - i)); i = j
        elif c == '/' and i + 1 < n and src[i + 1] == '*':
            j = src.find('*/', i + 2)
            if j < 0: j = n - 2
            seg = src[i:j + 2]
            out.append(re.sub(r'[^\n]', ' ', seg)); i = j + 2
        elif c == '"' or c == "'":
            q = c; j = i + 1
            while j < n and src[j] != q:
                if src[j] == '\\': j += 1
                if j < n and src[j] == '\n' and q == "'": break
                j += 1
            seg = src[i:j + 1]
            # keep the string content (it matters for io functions) but mark it
            out.append(seg); i = j + 1
        else:
            out.append(c); i += 1
    return ''.join(out)


TOK = re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d+\.?\d*(?:[eE][+-]?\d+)?|::|->|<<|>>|[<>=!+\-*/%&|^]=|&&|\|\||\+\+|--|\S')


def tokenize_norm(text):
    toks = []
    for t in TOK.findall(text):
        if t[0] in '"\'':
            toks.append('S')
        elif re.match(r'[A-Za-z_]', t):
            toks.append(t if t in KEYWORDS else 'ID')
        elif re.match(r'\d', t):
            toks.append('N')
        else:
            toks.append(t)
    return toks


def raw_tokens(text):
    return TOK.findall(text)


def find_functions(path):
    src = open(path, errors='replace').read()
    s = strip_comments(src)
    # strip preprocessor lines for brace matching? keep them, but they rarely contain braces.
    funcs = []
    n = len(s)
    depth = 0
    i = 0
    # stack of (type, startpos) where type in {'func','other'}
    stack = []
    line_of = [0] * (n + 1)
    ln = 1
    for k, ch in enumerate(s):
        line_of[k] = ln
        if ch == '\n': ln += 1
    line_of[n] = ln
    inside_func = 0
    while i < n:
        ch = s[i]
        if ch == '#' and (i == 0 or s[s.rfind('\n', 0, i) + 1:i].strip() == ''):
            j = s.find('\n', i)
            while j > 0 and s[j - 1] == '\\':
                j = s.find('\n', j + 1)
            i = n if j < 0 else j
            continue
        if ch == '{':
            kind = 'other'
            name = None
            if inside_func == 0:
                # look back for ')' optionally followed by const/noexcept/override/ctor init list
                head = None
                j = i - 1
                while j >= 0 and s[j].isspace(): j -= 1
                # skip trailing qualifiers / ctor init list: walk back to last ')' not past ; { }
                k = j
                # strip simple trailing words
                tail = s[max(0, j - 400):j + 1]
                mt = re.search(r'\)\s*(?:const|noexcept|override|final|\s)*$', tail)
                if not mt:
                    # ctor init list: find ') :' pattern
                    p = max(s.rfind(';', 0, i), s.rfind('}', 0, i), s.rfind('{', 0, i))
                    seg = s[p + 1:i]
                    mi = re.search(r'\)\s*(?:const)?\s*:[^:]', seg)
                    if mi and '::' not in seg[mi.start():mi.start()+3]:
                        k = p + 1 + mi.start()
                    else:
                        k = -1
                else:
                    k = j - len(tail) + 1 + mt.start()
                if k >= 0 and s[k] == ')':
                    # walk back to matching '('
                    d = 0; q = k
                    while q >= 0:
                        if s[q] == ')': d += 1
                        elif s[q] == '(':
                            d -= 1
                            if d == 0: break
                        elif s[q] in ';{}' : q = -1; break
                        q -= 1
                    if q > 0:
                        head = s[max(0, q - 300):q].rstrip()
                if head:
                    # name = last identifier (with qualifiers) of head
                    mm = re.search(r'((?:[\w~]+(?:<[^<>]*>)?::)*[~\w]+|operator\s*\S+)\s*$', head)
                    if mm:
                        cand = mm.group(1)
                        last = cand.split('::')[-1]
                        if last not in KEYWORDS and last not in ('if', 'for', 'while', 'switch', 'catch'):
                            kind = 'func'; name = cand; nameline = line_of[q]
                # in a class body (namespace/struct/class) we still allow functions; inside_func counts function nesting
            stack.append((kind, i, name, nameline if kind=='func' else 0))
            if kind == 'func':
                inside_func += 1
        elif ch == '}':
            if stack:
                kind, st, name, nameline = stack.pop()
                if kind == 'func':
                    inside_func -= 1
                    # find header start: go back to previous ';' or '}' or '{' at line granularity
                    hb = max(s.rfind(';', 0, st), s.rfind('}', 0, st), s.rfind('{', 0, st))
                    # skip preprocessor lines and blank lines
                    hstart = hb + 1
                    body = s[st:i + 1]
                    header = s[hstart:st]
                    # remove preprocessor lines from header
                    header = '\n'.join(l for l in header.split('\n') if not l.strip().startswith('#'))
                    funcs.append(dict(
                        file=os.path.relpath(path, ROOT), name=name,
                        start=nameline,
                        body_start=line_of[st], end=line_of[i],
                        header=' '.join(header.split()),
                        body=body,
                    ))
        i += 1
    return funcs


def main():
    allf = []
    for d in DIRS:
        for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
            for f in sorted(fn):
                if f.endswith(('.cpp', '.hpp', '.h')) and f != 'json.hpp' and 'Faddeeva' not in f:
                    allf.extend(find_functions(os.path.join(dp, f)))
    for f in allf:
        lines = [' '.join(l.split()) for l in f['body'].split('\n')]
        lines = [l for l in lines if l]
        f['lines'] = lines
        f['nlines'] = len(lines)
        f['ntoks'] = tokenize_norm(f['body'])
        f['rtoks'] = raw_tokens(f['body'])
        del f['body']
    json.dump(allf, open(OUT, 'w'))
    print(len(allf), 'functions;', sum(f['nlines'] for f in allf), 'body lines')


if __name__ == '__main__':
    main()
