from ipywidgets import HBox, VBox, Box,  Button, Layout, FloatProgress, Output
import hublib.ui as ui
import random, string
from string import Template

PARADIM = {}
PARADIM['SERVER'] = 'gateway2.marcc.jhu.edu'    
PARADIM['TUTORIAL_NAME'] = 'GaAs_Tutorial'
PARADIM['USER'] = 'dmejia4@jhu.edu'
PARADIM['CODE'] = ''
PARADIM['PWD'] = 'n4n0P4r4d1m.'
PARADIM['PORT'] = 22
PARADIM['SESSION'] = None
PARADIM['PBS_ID'] = 0
PARADIM['TIMER'] = 30
PARADIM['PW_TEMPLATE'] = Template('''&control
calculation = '${calculation}'
prefix = '${prefix}',
pseudo_dir = './',
outdir = './'
/
&system
ibrav = 2,
celldm(1) = 10.4749,
nat = 2,
ntyp = 2,
${system_params}
/
&electrons
/
ATOMIC_SPECIES
Ga 1.0 Ga.pz-bhs.UPF
As 1.0 As.pz-bhs.UPF
ATOMIC_POSITIONS ${crystal}
Ga 0.00 0.00 0.00
As 0.25 0.25 0.25
K_POINTS ${k_points}'''
)

PARADIM['EPSILON_TEMPLATE'] = Template('''&inputpp
outdir = './'
prefix = '${prefix}'
calculation = '${calculation}'
/
&energy_grid
smeartype = 'gauss'
intersmear = ${intersmear}
wmin = 0.0
wmax = 30.0
nw = 500
/'''
)

#SBATCH --reservation=$reservation
# this is not supported right now
PARADIM['PBS_TEMPLATE'] = Template('''#!/bin/bash -l
#SBATCH
#SBATCH --job-name=$job_name
#SBATCH --time=$time
#SBATCH --partition=$partition
#SBATCH --nodes=$nodes
#SBATCH --ntasks-per-node=$tasks
#SBATCH --mem-per-cpu=$mem
$modules
$command'''
)

def GetInputFile():
    in_input_file = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    return in_input_file

##################################################
# Conection
##################################################

PARADIM_UI = {}
PARADIM_UI['s0'] = {}
PARADIM_UI['s0']['pwd'] = ui.String(description='MARCC Password', name='MARCC Password',value=PARADIM['PWD'])
PARADIM_UI['s0']['code'] = ui.String(description='GOOGLE AUTH Code', name='GOOGLE AUTH Code',value=PARADIM['CODE'])
PARADIM_UI['s0']['user'] = ui.String(description='Marcc User name', name='Marcc User name',value=PARADIM['USER'])
PARADIM_UI['s0']['folder'] = ui.String(description='Working folder', name='Working folder',value=GetInputFile())
PARADIM_UI['s0']['button'] = Button(description='Connect')
PARADIM_UI['s0']['button'].w = Box([PARADIM_UI['s0']['button']])
PARADIM_UI['s0']['status'] = ui.String(name='Status', description='Status',value='')
PARADIM_UI['s0']['status'].dd.layout = Layout(width='100%')
PARADIM_UI['s0']['status'].disabled = True
PARADIM_UI['s0']['display'] = ui.Form([
            PARADIM_UI['s0']['user'], 
            PARADIM_UI['s0']['pwd'], 
            PARADIM_UI['s0']['code'], 
            PARADIM_UI['s0']['folder'], 
            PARADIM_UI['s0']['button'],
            PARADIM_UI['s0']['status'],            
        ], name = 'MARCC Credentials')
''',PARADIM_UI['s0']['commands']'''


##################################################
# First
##################################################

PARADIM_UI['s1'] = {}
PARADIM_UI['s1']['ecutwfc'] = ui.String(name='ecutwfc',description='ecutwfc',value='40.0')
PARADIM_UI['s1']['kpoints'] = ui.Text(
    description="kpoints",
    name="kpoints",
    value='''automatic
6 6 6 1 1 1'''    
)
PARADIM_UI['s1']['kpoints'].dd.layout = Layout(height='120px')
PARADIM_UI['s1']['input'] = ui.Text( description="inputdeck",  name="inputdeck", value='''''')
PARADIM_UI['s1']['input'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s1']['input'].disabled = True

PARADIM_UI['s1']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s1']['commands'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s1']['commands'].disabled = True

PARADIM_UI['s1']['stdin'] = ui.Text( description="stdin",name="stdin", value='''''')
PARADIM_UI['s1']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s1']['stdin'].disabled = True

PARADIM_UI['s1']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s1']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s1']['stdout'].disabled = True

PARADIM_UI['s1']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s1']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s1']['stderr'].disabled = True

PARADIM_UI['s1']['button'] = Button(description='Calculate Self Consistency')
PARADIM_UI['s1']['button'].layout = Layout(width='99%')
PARADIM_UI['s1']['button'].w = Box([PARADIM_UI['s1']['button']])

PARADIM_UI['s1']['job_id'] = ui.String(description='job ID',name='job ID',value='')
PARADIM_UI['s1']['job_id'].disabled = True
PARADIM_UI['s1']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s1']['status'].disabled = True
#PARADIM_UI['s1']['button_status'] = Button(description='Update Status')
#PARADIM_UI['s1']['button_status'].layout = Layout(width='99%')

PARADIM_UI['s1']['l1'] = VBox([
                               PARADIM_UI['s1']['ecutwfc'].w,
                               PARADIM_UI['s1']['kpoints'].w,
                               PARADIM_UI['s1']['button'].w,
                               PARADIM_UI['s1']['job_id'].w,
                               PARADIM_UI['s1']['status'].w,                               
                              ])
PARADIM_UI['s1']['l2'] = VBox([
                               PARADIM_UI['s1']['input'].w,
                               PARADIM_UI['s1']['commands'].w,                               
                              ])
PARADIM_UI['s1']['bs'] = HBox([PARADIM_UI['s1']['l1'],PARADIM_UI['s1']['l2']])
PARADIM_UI['s1']['bs'].w = Box([PARADIM_UI['s1']['bs']])
PARADIM_UI['s1']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s1']['bs'].layout = Layout(width='100%', border='1px')

s1_tab0 = ui.Form([PARADIM_UI['s1']['bs']], name = 'Crystal Inputs')
s1_tab1 = ui.Form([PARADIM_UI['s1']['stdin']], name = 'stdin')
s1_tab2 = ui.Form([PARADIM_UI['s1']['stdout']], name = 'stdout')
s1_tab3 = ui.Form([PARADIM_UI['s1']['stderr']], name = 'stderr')


def UpdateStep1( event ):
    global PARADIM_UI, PARADIM
    change = {  'k_points':PARADIM_UI['s1']['kpoints'].value, 
                'prefix':'gaas',
                'calculation':'scf',
                'crystal':'crystal',
                'system_params':'ecutwfc = ' + PARADIM_UI['s1']['ecutwfc'].value + "," }

    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s1']['input'].value = PARADIM['PW_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s1']['ecutwfc'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['kpoints'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['display'] = ui.Tab([s1_tab0, s1_tab1, s1_tab2, s1_tab3])
UpdateStep1({'type':'change','name':'value','new':'new'})

##################################################
# Second
##################################################

PARADIM_UI['s2'] = {}

PARADIM_UI['s2']['nbnd'] = ui.String(name='nbnd',description='nbnd',value='8')
PARADIM_UI['s2']['ecutwfc'] = ui.String(name='ecutwfc',description='ecutwfc',value='40.0')
PARADIM_UI['s2']['kpoints'] = ui.Text(
    description="kpoints",
    name="kpoints",
    value='''tpiba_b
3
0.500 0.500 0.500 20
0.000 0.000 0.000 20
1.000 0.000 0.000 20'''    
)
PARADIM_UI['s2']['kpoints'].dd.layout = Layout(height='120px')
PARADIM_UI['s2']['input'] = ui.Text( description="inputdeck",  name="inputdeck", value='''''')
PARADIM_UI['s2']['input'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s2']['input'].disabled = True

PARADIM_UI['s2']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s2']['commands'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s2']['commands'].disabled = True

PARADIM_UI['s2']['stdin'] = ui.Text( description="stdin", name="stdin", value='''''')
PARADIM_UI['s2']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stdin'].disabled = True

PARADIM_UI['s2']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s2']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stdout'].disabled = True

PARADIM_UI['s2']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s2']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stderr'].disabled = True

PARADIM_UI['s2']['button'] = Button(description='Calculate BandStructure')
PARADIM_UI['s2']['button'].layout = Layout(width='99%')
PARADIM_UI['s2']['button'].w = Box([PARADIM_UI['s2']['button']])

PARADIM_UI['s2']['job_id'] = ui.String(description='job ID',name='job ID',value='')
PARADIM_UI['s2']['job_id'].disabled = True
PARADIM_UI['s2']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s2']['status'].disabled = True

PARADIM_UI['s2']['l1'] = VBox([
                               PARADIM_UI['s2']['nbnd'].w,
                               PARADIM_UI['s2']['ecutwfc'].w,
                               PARADIM_UI['s2']['kpoints'].w,
                               PARADIM_UI['s2']['button'].w,
                               PARADIM_UI['s2']['job_id'].w,
                               PARADIM_UI['s2']['status'].w,                               
                              ])
PARADIM_UI['s2']['l2'] = VBox([
                               PARADIM_UI['s2']['input'].w,
                               PARADIM_UI['s2']['commands'].w,                               
                              ])
PARADIM_UI['s2']['bs'] = HBox([PARADIM_UI['s2']['l1'],PARADIM_UI['s2']['l2']])
PARADIM_UI['s2']['bs'].w = Box([PARADIM_UI['s2']['bs']])

PARADIM_UI['s2']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s2']['bs'].layout = Layout(width='100%', border='1px')

s2_tab0 = ui.Form([PARADIM_UI['s2']['bs']], name = 'band-structure Inputs')
s2_tab1 = ui.Form([PARADIM_UI['s2']['stdin']], name = 'stdin')
s2_tab2 = ui.Form([PARADIM_UI['s2']['stdout']], name = 'stdout')
s2_tab3 = ui.Form([PARADIM_UI['s2']['stderr']], name = 'stderr')


def UpdateStep2( event ):
    global PARADIM_UI, PARADIM
    change = {  'k_points':PARADIM_UI['s2']['kpoints'].value, 
                'prefix':'gaas',
                'calculation':'bands',
                'crystal':'',
                'system_params':'ecutwfc = ' + PARADIM_UI['s2']['ecutwfc'].value + ',\nnbnd = ' + PARADIM_UI['s2']['nbnd'].value + ',' 
             }

    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s2']['input'].value = PARADIM['PW_TEMPLATE'].safe_substitute(change)

PARADIM_UI['s2']['nbnd'].dd.observe(UpdateStep2)
PARADIM_UI['s2']['kpoints'].dd.observe(UpdateStep2)
UpdateStep2({'type':'change','name':'value','new':'new'})
PARADIM_UI['s2']['display'] = ui.Tab([s2_tab0, s2_tab1, s2_tab2, s2_tab3])


##################################################
# Third
##################################################

PARADIM_UI['s3'] = {}
PARADIM_UI['s3']['button'] = Button(description='Extract BandStructure')
PARADIM_UI['s3']['button'].layout = Layout(width='99%')
PARADIM_UI['s3']['button'].w = Box([PARADIM_UI['s3']['button']])

PARADIM_UI['s3']['input'] = ui.Text( description="inputs", name="inputs", value='''''')
PARADIM_UI['s3']['input'].dd.layout = Layout(width='100%', height='470px')
PARADIM_UI['s3']['input'].disabled = True

PARADIM_UI['s3']['output'] = Output()
PARADIM_UI['s3']['output'].layout = Layout(width='99%')

PARADIM_UI['s3']['l1'] = VBox([
                               PARADIM_UI['s3']['input'].w,
                               PARADIM_UI['s3']['button'].w,
                              ])
PARADIM_UI['s3']['l2'] = VBox([
                               PARADIM_UI['s3']['output'],
                              ])
PARADIM_UI['s3']['bs'] = HBox([PARADIM_UI['s3']['l1'],PARADIM_UI['s3']['l2']])
PARADIM_UI['s3']['bs'].w = Box([PARADIM_UI['s3']['bs']])

PARADIM_UI['s3']['l1'].layout = Layout(width='50%', border='1px')
PARADIM_UI['s3']['l2'].layout = Layout(width='50%', border='1px')
PARADIM_UI['s3']['bs'].layout = Layout(width='100%', border='1px')

PARADIM_UI['s3']['display'] = ui.Form([PARADIM_UI['s3']['bs']], name = 'band-structure Visualization')

##################################################
# Fourth
##################################################

PARADIM_UI['s4'] = {}

PARADIM_UI['s4']['nbnd'] = ui.String(description='nbnd',name='nbnd',value='16')
PARADIM_UI['s4']['ecutwfc'] = ui.String(description='ecutwfc',name='ecutwfc',value='40.0')
PARADIM_UI['s4']['intersmear'] = ui.String(description='intersmear',name='intersmear',value='0.2')
PARADIM_UI['s4']['kpoints'] = ui.Text(
    description="kpoints",
    name="kpoints",
    value='''automatic
10 10 10 1 1 1'''    
)
PARADIM_UI['s4']['kpoints'].dd.layout = Layout(height='120px')
PARADIM_UI['s4']['inputa'] = ui.Text( description="pw.x", name="pw.x", value='''''')
PARADIM_UI['s4']['inputa'].dd.layout = Layout(width='90%', height='180px')
PARADIM_UI['s4']['inputa'].disabled = True

PARADIM_UI['s4']['inputb'] = ui.Text( description="epsilon.x", name="epsilon.x", value='''''')
PARADIM_UI['s4']['inputb'].dd.layout = Layout(width='90%', height='180px')
PARADIM_UI['s4']['inputb'].disabled = True

PARADIM_UI['s4']['commands'] = ui.Text( description="commands", name="commands", value='''''')
PARADIM_UI['s4']['commands'].dd.layout = Layout(width='100%', height='150px')
PARADIM_UI['s4']['commands'].disabled = True

PARADIM_UI['s4']['stdin'] = ui.Text( description="stdin", name="stdin", value='''''')
PARADIM_UI['s4']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stdin'].disabled = True

PARADIM_UI['s4']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s4']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stdout'].disabled = True

PARADIM_UI['s4']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s4']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stderr'].disabled = True

PARADIM_UI['s4']['button'] = Button(description='Calculate Dielectric')
PARADIM_UI['s4']['button'].layout = Layout(width='99%')
PARADIM_UI['s4']['button'].w = Box([PARADIM_UI['s4']['button']])

PARADIM_UI['s4']['job_id'] = ui.String(description='job ID',name='job ID',value='')
PARADIM_UI['s4']['job_id'].disabled = True
PARADIM_UI['s4']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s4']['status'].disabled = True

PARADIM_UI['s4']['l1'] = VBox([
                               PARADIM_UI['s4']['nbnd'].w,
                               PARADIM_UI['s4']['ecutwfc'].w,
                               PARADIM_UI['s4']['intersmear'].w,
                               PARADIM_UI['s4']['kpoints'].w,
                               PARADIM_UI['s4']['button'].w,
                               PARADIM_UI['s4']['job_id'].w,
                               PARADIM_UI['s4']['status'].w,                               
                              ])
PARADIM_UI['s4']['l2'] = VBox([
                               PARADIM_UI['s4']['inputa'].w,
                               PARADIM_UI['s4']['inputb'].w,                               
                              ])
PARADIM_UI['s4']['bs'] = HBox([PARADIM_UI['s4']['l1'],PARADIM_UI['s4']['l2']])
PARADIM_UI['s4']['bs'].w = Box([PARADIM_UI['s4']['bs']])
PARADIM_UI['s4']['bs'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s4']['l2'].layout = Layout(width='100%', border='1px')

s4_tab0 = ui.Form([PARADIM_UI['s4']['bs']], name = 'band-structure Inputs')
s4_tab1 = ui.Form([PARADIM_UI['s4']['stdin']], name = 'stdin')
s4_tab2 = ui.Form([PARADIM_UI['s4']['stdout']], name = 'stdout')
s4_tab3 = ui.Form([PARADIM_UI['s4']['stderr']], name = 'stderr')


def UpdateStep4a( event ):
    global PARADIM_UI, PARADIM
    system_params = 'ecutwfc = ' + PARADIM_UI['s4']['ecutwfc'].value + ',\n'
    system_params = system_params + 'nbnd = ' + PARADIM_UI['s4']['nbnd'].value + ',\n'
    system_params = system_params + 'nosym = .true.\n'
    system_params = system_params + 'noinv = .true.'
    change = {  'k_points':PARADIM_UI['s4']['kpoints'].value, 
                'prefix':'gaas',
                'calculation':'nscf',
                'crystal':'',
                'system_params':system_params }
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s4']['inputa'].value = PARADIM['PW_TEMPLATE'].safe_substitute(change)

def UpdateStep4b( event ):
    global PARADIM_UI, PARADIM
    system_params = ''
    change = {  'intersmear':PARADIM_UI['s4']['intersmear'].value, 
                'prefix':'gaas',
                'calculation':'eps',
                'system_params':system_params }
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s4']['inputb'].value = PARADIM['EPSILON_TEMPLATE'].safe_substitute(change)
        
        
PARADIM_UI['s4']['nbnd'].dd.observe(UpdateStep4a)
PARADIM_UI['s4']['kpoints'].dd.observe(UpdateStep4a)
PARADIM_UI['s4']['intersmear'].dd.observe(UpdateStep4b)
UpdateStep4a({'type':'change','name':'value','new':'new'})
UpdateStep4b({'type':'change','name':'value','new':'new'})
PARADIM_UI['s4']['display'] = ui.Tab([s4_tab0, s4_tab1, s4_tab2, s4_tab3])



##################################################
# Fifth
##################################################


PARADIM_UI['s5'] = {}
PARADIM_UI['s5']['button'] = Button(description='Extract Coeficient')
PARADIM_UI['s5']['button'].layout = Layout(width='99%')
PARADIM_UI['s5']['button'].w = Box([PARADIM_UI['s5']['button']])

PARADIM_UI['s5']['input'] = ui.Text( description="inputs", name="inputs", value='''''')
PARADIM_UI['s5']['input'].dd.layout = Layout(width='100%', height='470px')
PARADIM_UI['s5']['input'].disabled = True

PARADIM_UI['s5']['output'] = Output()
PARADIM_UI['s5']['output'].layout = Layout(width='99%')

PARADIM_UI['s5']['l1'] = VBox([
                               PARADIM_UI['s5']['input'].w,
                               PARADIM_UI['s5']['button'].w,
                              ])
PARADIM_UI['s5']['l2'] = VBox([
                               PARADIM_UI['s5']['output'],
                              ])
PARADIM_UI['s5']['bs'] = HBox([PARADIM_UI['s5']['l1'],PARADIM_UI['s5']['l2']])
PARADIM_UI['s5']['bs'].w = Box([PARADIM_UI['s5']['bs']])

PARADIM_UI['s5']['l1'].layout = Layout(width='50%', border='1px')
PARADIM_UI['s5']['l2'].layout = Layout(width='50%', border='1px')
PARADIM_UI['s5']['bs'].layout = Layout(width='100%', border='1px')

PARADIM_UI['s5']['display'] = ui.Form([PARADIM_UI['s5']['bs']], name = 'band-structure Visualization')