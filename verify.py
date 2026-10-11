#!/usr/bin/env python3
"""Independent exact checker for the certificates in the Erdos-302 note.
Usage: python3 verify.py box_2_3.json [more.json ...]
Uses only the standard library (fractions, json, itertools)."""
import json, sys, itertools
from fractions import Fraction

def check(path):
    d=json.load(open(path)); P=d["primes"]; dims=d["dims"]; SC=d["scale"]
    cs=set()
    for e in itertools.product(*[range(x) for x in dims]):
        v=1
        for p,k in zip(P,e): v*=p**k
        cs.add(v)
    prev_th=0
    for i,j in enumerate(d["jumps"],1):
        assert j["k"]==i, "k not consecutive"
        th=j["threshold"]; assert th>=prev_th; prev_th=th
        load={}; tot=0
        for (a,b,c),w in j["packing"]:
            assert a in cs and b in cs and c in cs, "not a cell"
            assert len({a,b,c})==3, "not distinct"
            assert b*c==a*(b+c), "not harmonic: 1/a != 1/b+1/c"
            assert max(a,b,c)<=th, "triple above threshold"
            assert w>0
            for x in (a,b,c): load[x]=load.get(x,0)+w
            tot+=w
        assert all(v<=SC for v in load.values()), "vertex overloaded"
        assert Fraction(tot,SC)>j["k"]-1, "packing value too small"
    pref=Fraction(1)
    for p,a in zip(P,dims): pref*=Fraction(p-1,p)/(1-Fraction(1,p**a))
    D=pref*sum(Fraction(1,j["threshold"]) for j in d["jumps"])
    C=1-D
    print(path,"primes",P,"box",tuple(dims),"jumps",len(d["jumps"]),"OK")
    print("  deletion density D =",D)
    print("  f(N)/N <= %.12f + o(1)"%float(C))
    return d

if __name__=="__main__":
    for p in sys.argv[1:]: check(p)
