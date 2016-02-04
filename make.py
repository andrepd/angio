#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import count, product
import os

def frange(start, stop, step):
    for i in range(int((stop-start)/step)):
        yield start+i*step
    #r = start
    #while r < stop:
        #yield r
        #r += step

def permutations(d):
    '''Devolve uma lista com d,t para todas as combinações de parâmetros no dicionário d e para todas as listas de tips em t'''
    for i in (dict(zip(d, v)) for v in product(*d.values())):
        yield i
        #for j in t:
            #build(i,j)
            #run()
            #yield i,j

consts = {
        'DEecm': frange(0.2,1.2,0.05),
}

pwd = os.getcwd()
gcc = 'g++ main.cpp -std=c++11 -O3 -march=native -o main '
stuff = ['main.cpp', 'srj.hpp', 'run.sh']

# Criar pasta Run #N
N = 0
for i in count(1):
    if not os.path.exists('Run #{}'.format(i)):
        os.makedirs('Run #{}'.format(i))
        N = i
        break

# Mover stuff para lá
for i in stuff:
    os.system('cp {} Run\\ #{}'.format(i,N))
os.chdir('./Run #{}'.format(N))
print('Starting Run #{}'.format(N))

# Escrever cenas
with open('tasks.txt','w') as f:
    for i,d in zip(count(1), permutations(consts)):
        print('{}: '.format(i) + ', '.join('{} = {}'.format(a,b) for a,b in d.items()))

        #print(d)
        os.makedirs(str(i))
        comp = gcc + ' '.join('-{}={}'.format(a,b) for a,b in d.items())
        #print(comp)
        os.system(comp)
        
        with open('consts.txt', 'w') as g:
            g.write('\n'.join('{} = {}'.format(a,b) for a,b in d.items()))
        f.write('./{}/main >out.txt 2>err.txt\n'.format(i))
        
        os.system('mv main consts.txt {}'.format(i))
        os.system('cp ../plot.gp {}'.format(i))
        os.system('mkdir {}/data'.format(i))
