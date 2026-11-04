import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import sys,difflib
sys.path.insert(0, HERE)
from extract import strip_comments
def load(spec):
    f,rng=spec.split(':') if ':' in spec else (spec,None)
    L=strip_comments(open(f).read()).split('\n')
    if rng:
        a,b=map(int,rng.split('-')); L=L[a-1:b]
    return [' '.join(l.split()) for l in L if l.strip() and l.strip() not in ('{','}','};')]
a=load(sys.argv[1]); b=load(sys.argv[2])
sm=difflib.SequenceMatcher(None,a,b,autojunk=False)
common=sum(m.size for m in sm.get_matching_blocks())
print(f"{sys.argv[1]} ({len(a)}) vs {sys.argv[2]} ({len(b)}): common lines {common}, ratio {sm.ratio():.2f}")
