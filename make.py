#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import count
import os

def frange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

pre = '/home/angio/samm/angio/'
with open('tasks.txt','w') as f:
    for n,i in zip(count(), frange(0.2,1.2,0.05)):
        if not os.path.exists(str(i)):
            os.makedirs(str(i))
        comp = 'g++ main.cpp -std=c++11 -O3 -march=native -o main{0} -Doutput=\'\"{1}\"\' -DEecm={1}'.format(n,i)
        print comp
        os.system(comp)
        f.write('./main{}\n'.format(n))
