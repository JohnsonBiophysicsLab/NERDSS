import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import json, re, os, subprocess, collections

F=json.load(open(_os.path.join(OUTDIR, 'funcs.json')))
# gather all text (comment-stripped) per file
import sys; sys.path.insert(0,'.')
from extract import strip_comments
texts={}
for d in ['src','EXEs','include','run_code_tests']:
    for dp,dn,fn in os.walk(os.path.join(ROOT,d)):
        for f in fn:
            if f.endswith(('.cpp','.hpp','.h')) and f!='json.hpp':
                p=os.path.join(dp,f); texts[os.path.relpath(p,ROOT)]=strip_comments(open(p,errors='replace').read()).split('\n')
names=collections.defaultdict(list)
for f in F:
    last=f['name'].split('::')[-1]
    names[last].append(f)
for last,defs in sorted(names.items()):
    if last in ('main',) or last.startswith('operator') or '~' in last: continue
    pat=re.compile(r'(?<![\w])'+re.escape(last)+r'\s*(\(|<)')
    pat2=re.compile(r'(?<![\w.>:])&?'+re.escape(last)+r'(?![\w])')
    uses=[]
    deflines={(d['file'],d['start']) for d in defs}
    for fn,lines in texts.items():
        for i,l in enumerate(lines,1):
            if (fn,i) in deflines: continue
            if re.search(r'(?<![\w])'+re.escape(last)+r'(?![\w])', l):
                # skip declarations: line in header ending with ; that has a return type before name
                s=l.strip()
                isdecl = re.match(r'^(static\s+|inline\s+|virtual\s+|extern\s+|friend\s+)*[\w:<>,\s\*&]+[\s\*&]'+re.escape(last)+r'\s*\(', s) and not re.search(r'=\s*'+re.escape(last), s) and not s.startswith('return')
                if isdecl: continue
                uses.append((fn,i,s[:100]))
    if not uses:
        for d in defs:
            print(f"DEAD? {d['file']}:{d['start']} {d['name']} ({d['nlines']} lines)")
