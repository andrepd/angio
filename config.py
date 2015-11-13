#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Ferramenta para configurar, compilar e correr o programa. 
  Se o script for corrido, lê as definições do ficheiro, usa-as para compilar o programa, e corre
  Se o script for importado, define as funções parse, config e run
'''
import os
from sys import argv, exit
from time import clock, time
from itertools import product

defaults = {
        'Lx': '128',
        'Ly': '128',
        'medini': '-.2',
        'Amp': '10.55',
        'D': '11',
        'dt': '.02',
        'rad': '5',
        'vbase': '.01',
        'valoralfa': '.065',
        'nmax': '4',
        'tf': '100000',
        'vmax':'0.15',
	'Pmax':'0.03',
	'M1':'1',
	'rho0':'1', 
	'Eecm': '0.6329',
	'Ec':'0.084388',
	'nuc':'0.49',
	'nuecm':'0.13',
	'vconc':'1',
	'output':'data',
        'passo': '100',
        'passotips': '500',
}

def box(s,c='*',n=None):
    if n == None:
        n = max([len(x) for x in s])+6
    return c*n+'\n' + '\n'.join([c+x.center(n-2)+c for x in s]) + '\n'+c*n+'\n'

def parse(name):
    '''Lê as definições do ficheiro name e retorna um dicionário com os defines e uma lista com as coordenadas iniciais das tip cells'''
    print 'Reading input file...',
    defines = {}
    tips = []
    with open(name, 'r') as f:
        f = [x for x in f if x[0] != '#' and x != '\n' and x != 'TIPS\n']
        for i in f:
            if i[0] == '\t':
                t = i[1:].split()
                tips.append((int(i[0]), int(i[1])))
            else:
                x = i.strip(' ').strip('\t').strip('\n').split('=')
                defines[x[0]] = x[1]
    print 'Done.'
    return defines,tips

def build(defines, tips):
    '''Usando os defines em defines e as coordenadas das tips em tips compila o programa'''
    for i in defaults:
        if i not in defines:
            defines[i] = defaults[i]

    if not tips:
        tips.append((int(defines['Lx'])/5+10,int(defines['Ly'])/2))
    with open('tips.in','w') as f:
        for a,b in tips:
            f.write(str(a)+' '+str(b)+'\n')

    if not os.path.exists(defines['output']):
        os.makedirs(defines['output'])
    defines['output'] = '\'\"'+defines['output']+'\"\''

    print 'Constants:'
    for i in defines:
        print i,'=',defines[i]
    print

    print 'Tips:'
    for i in tips:
        print i
    print

    print 'Compiling...'
    #os.system('ulimit -s '+str(int(defines['Lx'])*int(defines['Ly'])*128/1024))
    os.system('g++ -pg main.cpp -std=c++11 -march=native -O3 -o main '+' '.join(['-D'+i+'='+j for i,j in defines.items()]))
    print 'Done.\n'

def run():
    '''Corre o programa, contabilizando o tempo'''
    print 'Running:\n'
    print '-----\n'
    ti = time()
    os.system('./main')
    tf = time()
    print '\n-----\n'
    print 'Running time:',round(tf-ti,2),'seconds.'
    print 'Done.'

def permutations(d,t):
    '''Devolve uma lista com d,t para todas as combinações de parâmetros no dicionário d e para todas as listas de tips em t'''
    for i in [dict(zip(d, v)) for v in product(*d.values())]:
        for j in t:
            #build(i,j)
            #run()
            yield i,j

if __name__ == '__main__':
    if len(argv) != 2:
        print 'Usage:',argv[0],'file'
        exit()

    print box([
    'Computational Biology',
    'University of Coimbra - Portugal',
    ' ',
    'Sprouting Angiogenesis',
    ' ',
    'Phase-field model with elasticity',
    ' ',
    'Developers:',
    'Rui Travasso, PhD',
    'Andre Duarte, BSc',
    'Marcos Gouveia, BSc',
    'Marina Oliveira, Student'
    ])

    d,t = parse(argv[1])
    build(d,t)
    run()

