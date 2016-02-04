#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import count
from os import *

consts = {
        'DEecm': frange(0.2,1.2,0.05),
}

def frange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def permutations(d,t):
    '''Devolve uma lista com d,t para todas as combinações de parâmetros no dicionário d e para todas as listas de tips em t'''
    for i in (zip(d, v) for v in product(*d.values())):
        yield i
        #for j in t:
            #build(i,j)
            #run()
            #yield i,j

pwd = getcwd()
gcc = 'g++ main.cpp -std=c++11 -O3 - march=native -o main '
stuff = ['main.cpp', 'srj.hpp', 'run.sh']

# Criar Run #N e mover para lá stuff

N = 0
for i in count(1):
    if not path.exists('Run #{}'.format(i)):
        makedirs('Run #{}'.format(i))
        N = i
        break

for i in stuff:
    system('cp {} Run\\ #{}'.format(i,N))
chdir('./Run #{}'.format(N))
print('Starting Run #{}'.format(N))

with open('tasks.txt','w') as f:
    for i,d in zip(count(1), permutations(consts)):
        makedirs(str(i))
        comp = gcc + ['-{}={}'.format(a,b) for a,b, in consts].join(' ')
        print comp
        #system(comp)
        f.write('./{}main\n'.format(i))
        system('mv main {}'.format(i))
