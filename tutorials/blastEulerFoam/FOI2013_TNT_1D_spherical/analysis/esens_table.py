import numpy as np, glob, re, sys
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz
p0=101325.0
def probe(c):
    fs=sorted(glob.glob(f"{c}/postProcessing/pressureProbes/*/p"),key=lambda f:float(f.split('/')[-2]))
    xs=[float(re.search(r'\(([-\d.e]+)',l).group(1)) for l in open(fs[0]) if l.startswith('# Probe')]
    d=np.vstack([np.loadtxt(f,comments='#',ndmin=2) for f in fs]); d=d[np.argsort(d[:,0],kind='stable')]
    _,u=np.unique(d[:,0],return_index=True); return xs,d[u]
res={}
for c in sys.argv[1:]:
    xs,d=probe(c); t=d[:,0]
    for r in (2.828,4.472,6.325):
        i=xs.index(r); dp=d[:,1+i]-p0
        if dp.max()<0.05*p0: res[c,r]=None; continue
        k=np.argmax(dp); ia=np.argmax(dp>0.05*dp[k]); neg=np.where(dp[k:]<0)[0]
        ie=k+neg[0] if len(neg) else len(t)
        res[c,r]=(t[ia]*1e3,dp[k]/1e3,np.trapezoid(dp[ia:ie],t[ia:ie]),len(neg)>0)
for r in (2.828,4.472,6.325):
    print(f"r = {r} m")
    b=res[sys.argv[1],r]
    for c in sys.argv[1:]:
        v=res[c,r]
        if v is None: print(f"  {c}: not reached"); continue
        s=f"  {c}: ta {v[0]:.3f} ms  pk {v[1]:6.1f} kPa  I+ {v[2]:6.1f} Pa s" + ("" if v[3] else " (open)")
        if b and c!=sys.argv[1] and b[3] and v[3]: s+=f"   | dta {v[0]-b[0]:+.3f} ms, pk {100*(v[1]/b[1]-1):+.1f}%, I+ {100*(v[2]/b[2]-1):+.1f}%"
        print(s)
