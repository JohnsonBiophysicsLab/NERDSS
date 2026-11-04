import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import json, os, re, collections, sys
sys.path.insert(0,'.')
from extract import find_functions

F=json.load(open(_os.path.join(OUTDIR, 'funcs.json')))
# also add run_code_tests functions as roots
roots_extra=[]
for fn in os.listdir(ROOT+'/run_code_tests'):
    if fn.endswith('.cpp'):
        for f in find_functions(os.path.join(ROOT,'run_code_tests',fn)):
            roots_extra.append(set(re.findall(r'[A-Za-z_]\w*', f['body'])))
byname=collections.defaultdict(list)
for i,f in enumerate(F): byname[f['name'].split('::')[-1]].append(i)
def mentions(i): return set(t for t in F[i]['rtoks'] if re.match(r'[A-Za-z_]',t))
M=[mentions(i) for i in range(len(F))]
seen=set(); stack=[]
for i,f in enumerate(F):
    if f['file'].startswith('mpi_proto'): continue
    if f['name']=='main' or '::' in f['name'] or f['name'].startswith('operator'):
        stack.append(i)
# global-scope initializers in non-function code: crude - add names mentioned outside functions? skip
for s in roots_extra:
    for n in s:
        stack.extend(byname.get(n,[]))
while stack:
    i=stack.pop()
    if i in seen: continue
    seen.add(i)
    for n in M[i]:
        for j in byname.get(n,[]):
            if j not in seen and not F[j]['file'].startswith('mpi_proto'): stack.append(j)
dead=[f for i,f in enumerate(F) if i not in seen and not f['file'].startswith('mpi_proto')]
for f in sorted(dead,key=lambda f:(f['file'],f['start'])):
    print(f"{f['file']}:{f['start']} {f['name']} ({f['nlines']} lines)")
print(len(dead),'unreachable;',sum(f['nlines'] for f in dead),'lines')
