#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Exemplo de script para farmar'''

from config import *

consts = {
    'valoralfa': ['.060','.065','.070'],
    'passotips': [str(x) for x in range(100,1000,200)]
}

tips = [
    [],
    [(1,2),(2,4)]
]

for i,j in permutations(consts, tips):
    print i,j
