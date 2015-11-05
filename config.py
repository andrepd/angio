import os
from sys import argv, exit
from time import clock, time

defaults = {
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
        'vmax':'0.15',
	'Pmax':'0.03',
	'L0':'1',
	'M1':'1',
	'rho0':'1',
	'vconc':'1',
	'output':'\'\"data\"\''
}

def box(s,c='*',n=None):
    if n == None:
        n = max([len(x) for x in s])+6
    return c*n+'\n' + '\n'.join([c+x.center(n-2)+c for x in s]) + '\n'+c*n+'\n'

def parse(name):
    print 'Reading input file...',
    defines = defaults
    tips = []
    with open(name, 'r') as f:
        f = [x for x in f if x[0] != '#' and x != '\n' and x != 'TIPS\n']
        for i in f:
            if i[0] == '\t':
                tips.append(i[1:])
            else:
                x = i.strip(' ').strip('\t').strip('\n').split('=')
                defines[x[0]] = x[1]
    print 'Done.'
    return defines,tips

def build(defines, tips):
    if not tips:
        tips.append(str(int(defines['Lx'])/5+10)+' '+str(int(defines['Ly'])/2))
    with open('tips.in','w') as f:
        for i in tips:
            f.write(i)

    d['output'] = '\'\"'+d['output']+'\"\''

    print 'Constants:'
    for i in defines:
        print i,'=',defines[i]
    print

    print 'Tips:'
    for i in tips:
        print i
    print

    if not os.path.exists(d['output'][2:-2]):
        os.makedirs(d['output'][2:-2])

    print 'Compiling...'
    os.system('ulimit -s '+str(int(defines['Lx'])*int(defines['Ly'])*128/1024))
    os.system('g++ -pg main.cpp -std=c++11 -march=native -O3 -o main '+' '.join(['-D'+i+'='+j for i,j in defines.items()]))
    print 'Done.\n'

def run():
    print 'Running:\n'
    print '-----\n'
    ti = time()
    os.system('./main')
    tf = time()
    print '\n-----\n'
    print 'Running time:',round(tf-ti,2),'seconds.'
    print 'Done.'

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



