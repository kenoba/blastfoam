import numpy as np, re, sys, os
def F(t,n):
    s=open(f"{t}/{n}").read(); m=re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(",s)
    if m: return np.array(s[m.end():s.index(")",m.end())].split(),float)
    return float(re.search(r"internalField\s+uniform\s+([^;]+);",s).group(1))
def peak_profile(t, xs, ys):
    X,Y,P=F(t,"Ccx"),F(t,"Ccy"),F(t,"pTemporalMax")
    for x in xs:
        row=[]
        for y in ys:
            k=np.argmin((X-x)**2+(Y-y)**2); row.append((P[k]-101325)/1e3)
        print(f"  x={x}: "+" ".join(f"{v:6.1f}" for v in row))
t=sys.argv[1]
ys=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2.0]
print(f"peak overpressure [kPa] (pTemporalMax at {t}) at y =", ys)
peak_profile(t,[2.0,3.0,4.0,5.0,6.0],ys)
