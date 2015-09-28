from os import system
from sys import argv, exit

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

defines = []
tips = []
with open(argv[1], 'r') as f:
    f = [x for x in f if x[0] != '#' and x != '\n' and x != 'TIPS\n']
    for i in f:
        if i[0] == '\t':
            tips.append(i[1:])
        else:
            defines.append('-D '+i)

with open('tips.in','w') as f:
    for i in tips:
        f.write(i)

print defines
print tips

print 'Compiling... ',
system('g++ main.cpp -std=c++11 -march=native -O3 -o main'+' '.join(defines))
print 'Done.'

print 'Running:\n'
system('./main')
