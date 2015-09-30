from os import system
from sys import argv, exit
from time import clock,time

defines = {
        'Lx': '128',
        'Ly': '128',
        'medini': '-.2',
        'ge': '.126',
        'Amp': '10.55',
        'D': '11',
        'dt': '.02',
        'rad': '5',
        'vbase': '.01',
        'valoralfa': '.065',
        'nmax': '4',
        'iNf': '1',
        'tf': '100000',
}

if len(argv) != 2:
    print 'Usage:',argv[0],'file'
    exit()

def box(s,c='*',n=None):
    if n == None:
        n = max([len(x) for x in s])+6
    return c*n+'\n' + '\n'.join([c+x.center(n-2)+c for x in s]) + '\n'+c*n

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

print 'Reading input file...\n'

#defines = []
tips = []
with open(argv[1], 'r') as f:
    f = [x for x in f if x[0] != '#' and x != '\n' and x != 'TIPS\n']
    for i in f:
        if i[0] == '\t':
            tips.append(i[1:])
        else:
            x = i.strip(' ').strip('\t').strip('\n').split('=')
            defines[x[0]] = x[1]

if not tips:
    tips.append(str(int(defines['Lx'])/5+10)+' '+str(int(defines['Ly'])/2))
with open('tips.in','w') as f:
    for i in tips:
        f.write(i)

print 'Compiling... ',
system('ulimit -s '+str(int(defines['Lx'])*int(defines['Ly'])*128/1024))
system('g++ -pg main.cpp -std=c++11 -march=native -O3 -o main '+' '.join(['-D'+i+'='+j for i,j in defines.items()]))
print 'Done.'

print 'Running:\n'
print '-----\n'
ti = time()
system('./main')
tf = time()
print '\n-----\n'
print 'Running time:',round(tf-ti,2),'seconds.'
print 'Done.'
