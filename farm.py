#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Exemplo de script para farmar'''

from config import *
import os
import json

consts = {
    'valoralfa': ['.060','.065','.070'],
}

tips = [
    [],
]

n=1
for i,j in permutations(consts, tips):
    #print i,j
    build(dict(i),list(j))
    N = str(n)
    if not os.path.exists(N):
        os.makedirs(N)
    os.rename('main', N+'/main')
    os.rename('tips.in', N+'/tips.in')
    if not os.path.exists(N+'/data'):
        os.makedirs(N+'/data')
    with open(N+'/consts.txt','w') as f:
        f.write(json.dumps(i))
        f.write('\n')
        f.write(json.dumps(j))
    n+=1

with open('run.sh','w') as f:
    for i in range(1,n):
        I = str(i)
        f.write('cd '+I+' ; ./main\n')
