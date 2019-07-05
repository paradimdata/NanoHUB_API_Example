from ipywidgets import HBox, VBox, Box,  Button, Layout, FloatProgress, Output
import hublib.ui as ui
import random, string
from string import Template

PARADIM = {}
PARADIM['SERVER'] = 'gateway2.marcc.jhu.edu'    
PARADIM['TUTORIAL_NAME'] = 'LaNiO3_Tutorial'
PARADIM['USER'] = ''
PARADIM['CODE'] = ''
PARADIM['PWD'] = ''
PARADIM['PORT'] = 22
PARADIM['SESSION'] = None
PARADIM['PBS_ID'] = 0
PARADIM['TIMER'] = 30
PARADIM['PW_TEMPLATE'] = Template('''&control
 calculation = 'scf',
 prefix = '${prefix}',
 pseudo_dir = './',
 outdir = './',
 nstep = 300, 
/
&system
 ibrav = 0,
 celldm(1) = ${celldm},
 nat = ${nat},
 ntyp = ${ntyp},
 ecutwfc = ${ecutwfc},
 occupations = 'smearing',
 smearing    = 'gauss',
 degauss     = 1.00000e-02,
/
&electrons
 conv_thr = 1.0d-8
/
&IONS
/
&CELL
/

ATOMIC_SPECIES
  La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
  Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
  O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
CELL_PARAMETERS (alat=  7.28886267)
   0.973022005   0.000000000   0.000000000
   0.000000000   0.973022005   0.000000000
   0.000000000   0.000000000   0.973022005
  
ATOMIC_POSITIONS (crystal)
La       0.500000000   0.500000000   0.500000000
Ni       0.000000000   0.000000000   0.000000000
O        0.000000000   0.000000000   0.500000000
O        0.500000000   0.000000000   0.000000000
O        0.000000000   0.500000000   0.000000000
K_POINTS automatic
4 4 4 1 1 1
'''
)


PARADIM['PH_TEMPLATE'] = Template('''phonons of LaNiO3
&inputph
 prefix = '${prefix}',
 amass(1) = ${amass1},
 amass(2) = ${amass2},
 amass(3) = ${amass3},
 fildyn = ${fildyn},
 tr2_ph = 1.0d-14,
 ldisp = .true.,
 nq1 = ${nq1},
 nq2 = ${nq2},
 nq3 = ${nq3},
/
'''
)

PARADIM['Q2R_TEMPLATE'] = Template('''&input
 fildyn = '${fildyn}',
 flfrc = '${flfrc}'
/
'''
)

PARADIM['MATDYN_TEMPLATE'] = Template('''&input
 asr = 'simple',
 flfrc = '${flfrc}'
 flfrq = '${flfrq}'
/
${qpoints}
'''
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
PARADIM_UI['s1']['ecutwfc'] = ui.String(name='ecutwfc',description='ecutwfc',value='100')
PARADIM_UI['s1']['ntyp'] = ui.String(name='ntyp',description='ntyp',value='3')
PARADIM_UI['s1']['nat'] = ui.String(name='nat',description='nat',value='5')
PARADIM_UI['s1']['celldm'] = ui.String(name='celldm',description='celldm',value='7.28886267')

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

PARADIM_UI['s1']['button'] = Button(description='Calculate Diamond Structure')
PARADIM_UI['s1']['button'].layout = Layout(width='99%')
PARADIM_UI['s1']['button'].w = Box([PARADIM_UI['s1']['button']])

PARADIM_UI['s1']['job_id'] = ui.String(description='job ID',name='job ID',value='')
PARADIM_UI['s1']['job_id'].disabled = True
PARADIM_UI['s1']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s1']['status'].disabled = True


PARADIM_UI['s1']['l1'] = VBox([
                               PARADIM_UI['s1']['ecutwfc'],
                               PARADIM_UI['s1']['ntyp'],
                               PARADIM_UI['s1']['nat'],
                               PARADIM_UI['s1']['celldm'],
                               PARADIM_UI['s1']['button'].w,
                               PARADIM_UI['s1']['job_id'],
                               PARADIM_UI['s1']['status'],                               
                              ])
PARADIM_UI['s1']['l2'] = VBox([
                               PARADIM_UI['s1']['input'],
                               PARADIM_UI['s1']['commands'],
                              ])
PARADIM_UI['s1']['bs'] = HBox([PARADIM_UI['s1']['l1'],PARADIM_UI['s1']['l2']])
PARADIM_UI['s1']['bs'].w = Box([PARADIM_UI['s1']['bs']])
PARADIM_UI['s1']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s1']['l1'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s1']['bs'].layout = Layout(width='100%', border='1px')

s1_tab0 = ui.Form([PARADIM_UI['s1']['bs']], name = 'Crystal Inputs')
s1_tab1 = ui.Form([PARADIM_UI['s1']['stdin']], name = 'stdin')
s1_tab2 = ui.Form([PARADIM_UI['s1']['stdout']], name = 'stdout')
s1_tab3 = ui.Form([PARADIM_UI['s1']['stderr']], name = 'stderr')


def UpdateStep1( event ):
    global PARADIM_UI, PARADIM
    change = {  'k_points':'automatic \n 4 4 4 1 1 1\n', 
                'prefix':'LaNiO3',
                'calculation':'scf',
                'nat':PARADIM_UI['s1']['nat'].value,
                'ntyp':PARADIM_UI['s1']['ntyp'].value,
                'celldm':PARADIM_UI['s1']['celldm'].value,
                'ecutwfc':PARADIM_UI['s1']['ecutwfc'].value,
                'crystal':'crystal', 
                'cl_position': '0.25', 
    }

    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s1']['input'].value = PARADIM['PW_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s1']['ecutwfc'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['ntyp'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['nat'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['celldm'].dd.observe(UpdateStep1)
PARADIM_UI['s1']['display'] = ui.Tab([s1_tab0, s1_tab1, s1_tab2, s1_tab3])
UpdateStep1({'type':'change','name':'value','new':'new'})

##################################################
# Second
##################################################

PARADIM_UI['s2'] = {}
PARADIM_UI['s2']['amass1'] = ui.String(name='amass',description='amass',value='138.90550')
PARADIM_UI['s2']['amass2'] = ui.String(name='amass',description='amass',value='58.69340')
PARADIM_UI['s2']['amass3'] = ui.String(name='amass',description='amass',value='5.99940')
PARADIM_UI['s2']['nq'] = ui.String(name='nq(1-3)',description='nq(1-3)',value='2')

PARADIM_UI['s2']['input'] = ui.Text( description="inputdeck",  name="inputdeck", value='''''')
PARADIM_UI['s2']['input'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s2']['input'].disabled = True

PARADIM_UI['s2']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s2']['commands'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s2']['commands'].disabled = True

PARADIM_UI['s2']['stdin'] = ui.Text( description="stdin",name="stdin", value='''''')
PARADIM_UI['s2']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stdin'].disabled = True

PARADIM_UI['s2']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s2']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stdout'].disabled = True

PARADIM_UI['s2']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s2']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s2']['stderr'].disabled = True

PARADIM_UI['s2']['button'] = Button(description='Calculate Potential Frequencies')
PARADIM_UI['s2']['button'].layout = Layout(width='99%')
PARADIM_UI['s2']['button'].w = Box([PARADIM_UI['s2']['button']])

PARADIM_UI['s2']['job_id'] = ui.String(description='job ID',name='job ID',value='')
PARADIM_UI['s2']['job_id'].disabled = True
PARADIM_UI['s2']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s2']['status'].disabled = True


PARADIM_UI['s2']['l1'] = VBox([
                               PARADIM_UI['s2']['amass1'],
                               PARADIM_UI['s2']['amass2'],
                               PARADIM_UI['s2']['amass3'],
                               PARADIM_UI['s2']['nq'],
                               PARADIM_UI['s2']['button'].w,
                               PARADIM_UI['s2']['job_id'],
                               PARADIM_UI['s2']['status'],                               
                              ])
PARADIM_UI['s2']['l2'] = VBox([
                               PARADIM_UI['s2']['input'],
                               PARADIM_UI['s2']['commands'],                               
                              ])
PARADIM_UI['s2']['bs'] = HBox([PARADIM_UI['s2']['l1'],PARADIM_UI['s2']['l2']])
PARADIM_UI['s2']['bs'].w = Box([PARADIM_UI['s2']['bs']])
PARADIM_UI['s2']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s2']['l1'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s2']['bs'].layout = Layout(width='100%', border='1px')

s2_tab0 = ui.Form([PARADIM_UI['s2']['bs']], name = 'Crystal Inputs')
s2_tab1 = ui.Form([PARADIM_UI['s2']['stdin']], name = 'stdin')
s2_tab2 = ui.Form([PARADIM_UI['s2']['stdout']], name = 'stdout')
s2_tab3 = ui.Form([PARADIM_UI['s2']['stderr']], name = 'stderr')


def UpdateStep2( event ):
    global PARADIM_UI, PARADIM
    change = {  
                'prefix':'LaNiO3',
                'amass1' : PARADIM_UI['s2']['amass1'].value, 
                'amass2' : PARADIM_UI['s2']['amass2'].value, 
                'amass3' : PARADIM_UI['s2']['amass3'].value, 
                'fildyn': 'dyn',
                'nq1': PARADIM_UI['s2']['nq'].value,
                'nq2': PARADIM_UI['s2']['nq'].value,
                'nq3': PARADIM_UI['s2']['nq'].value
    }    
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s2']['input'].value = PARADIM['PH_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s2']['amass1'].dd.observe(UpdateStep2)
PARADIM_UI['s2']['amass2'].dd.observe(UpdateStep2)
PARADIM_UI['s2']['amass3'].dd.observe(UpdateStep2)
PARADIM_UI['s2']['nq'].dd.observe(UpdateStep2)
PARADIM_UI['s2']['display'] = ui.Tab([s2_tab0, s2_tab1, s2_tab2, s2_tab3])
UpdateStep2({'type':'change','name':'value','new':'new'})


##################################################
# Third
##################################################

PARADIM_UI['s3'] = {}
PARADIM_UI['s3']['flfrc'] = ui.String(name='flfrc',description='flfrc',value='diam.fc')

PARADIM_UI['s3']['input'] = ui.Text( description="inputdeck",  name="inputdeck", value='''''')
PARADIM_UI['s3']['input'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s3']['input'].disabled = True

PARADIM_UI['s3']['output'] = ui.Text( description="output",  name="output", value='''''')
PARADIM_UI['s3']['output'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s3']['output'].disabled = True

PARADIM_UI['s3']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s3']['commands'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s3']['commands'].disabled = True

PARADIM_UI['s3']['stdin'] = ui.Text( description="stdin",name="stdin", value='''''')
PARADIM_UI['s3']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s3']['stdin'].disabled = True

PARADIM_UI['s3']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s3']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s3']['stdout'].disabled = True

PARADIM_UI['s3']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s3']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s3']['stderr'].disabled = True

PARADIM_UI['s3']['button'] = Button(description='Calculate Interatomic Force')
PARADIM_UI['s3']['button'].layout = Layout(width='99%')
PARADIM_UI['s3']['button'].w = Box([PARADIM_UI['s3']['button']])

PARADIM_UI['s3']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s3']['status'].disabled = True

PARADIM_UI['s3']['l1'] = VBox([
                               PARADIM_UI['s3']['flfrc'].dd,
                               PARADIM_UI['s3']['button'].w,
                               PARADIM_UI['s3']['status'].dd
                              
                              ])
PARADIM_UI['s3']['l2'] = VBox([
                               PARADIM_UI['s3']['input'].dd,
                               PARADIM_UI['s3']['output'].dd,
                              ])
PARADIM_UI['s3']['bs'] = HBox([PARADIM_UI['s3']['l1'],PARADIM_UI['s3']['l2']])
PARADIM_UI['s3']['bs'].w = Box([PARADIM_UI['s3']['bs']])
PARADIM_UI['s3']['l1'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s3']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s3']['bs'].layout = Layout(width='100%', border='1px')

s3_tab0 = ui.Form([PARADIM_UI['s3']['bs']], name = 'Crystal Inputs')
s3_tab1 = ui.Form([PARADIM_UI['s3']['stdin']], name = 'stdin')
s3_tab2 = ui.Form([PARADIM_UI['s3']['stdout']], name = 'stdout')
s3_tab3 = ui.Form([PARADIM_UI['s3']['stderr']], name = 'stderr')


def UpdateStep3( event ):
    global PARADIM_UI, PARADIM
    change = {  
                'fildyn' : 'dyn',
                'flfrc' : PARADIM_UI['s3']['flfrc'].value,
             }
  
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s3']['input'].value = PARADIM['Q2R_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s3']['flfrc'].dd.observe(UpdateStep3)

PARADIM_UI['s3']['display'] = ui.Tab([s3_tab0, s3_tab1, s3_tab2, s3_tab3])
UpdateStep3({'type':'change','name':'value','new':'new'})



##################################################
# Fourth
##################################################

PARADIM_UI['s4'] = {}
PARADIM_UI['s4']['flfrq'] = ui.String(name='flfrq',description='flfrq',value='diam.freq')
PARADIM_UI['s4']['qpoints'] = ui.Text(
    description="qpoints",
    name="qpoints",
    value='''31
0.5000  0       0
0.5000  0.0500  0.0500
0.5000  0.1000  0.1000
0.5000  0.1500  0.1500
0.5000  0.2000  0.2000
0.5000  0.2500  0.2500
0.5000  0.3000  0.3000
0.5000  0.3500  0.3500
0.5000  0.4000  0.4000
0.5000  0.4500  0.4500
0.5000  0.5000  0.5000
0.5000  0.5000  0.4500
0.5000  0.5000  0.4000
0.5000  0.5000  0.3500
0.5000  0.5000  0.3000
0.5000  0.5000  0.2500
0.5000  0.5000  0.2000
0.5000  0.5000  0.1500
0.5000  0.5000  0.1000
0.5000  0.5000  0.0500
0.5000  0.5000  0
0.4500  0.4500  0
0.4000  0.4000  0
0.3500  0.3500  0
0.3000  0.3000  0
0.2500  0.2500  0
0.2000  0.2000  0
0.1500  0.1500  0
0.1000  0.1000  0
0.0500  0.0500  0
0       0       0
'''    
)
PARADIM_UI['s4']['qpoints'].dd.layout = Layout(height='120px')

PARADIM_UI['s4']['input'] = ui.Text( description="inputdeck",  name="inputdeck", value='''''')
PARADIM_UI['s4']['input'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s4']['input'].disabled = True

PARADIM_UI['s4']['output'] = ui.Text( description="output",  name="output", value='''''')
PARADIM_UI['s4']['output'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s4']['output'].disabled = True

PARADIM_UI['s4']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s4']['commands'].dd.layout = Layout(width='90%', height='150px')
PARADIM_UI['s4']['commands'].disabled = True

PARADIM_UI['s4']['stdin'] = ui.Text( description="stdin",name="stdin", value='''''')
PARADIM_UI['s4']['stdin'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stdin'].disabled = True

PARADIM_UI['s4']['stdout'] = ui.Text( description="stdout", name="stdout", value='''''')
PARADIM_UI['s4']['stdout'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stdout'].disabled = True

PARADIM_UI['s4']['stderr'] = ui.Text( description="stderr", name="stderr", value='''''')
PARADIM_UI['s4']['stderr'].dd.layout = Layout(width='100%', height='300px')
PARADIM_UI['s4']['stderr'].disabled = True

PARADIM_UI['s4']['button'] = Button(description='Calculate Interatomic Force')
PARADIM_UI['s4']['button'].layout = Layout(width='99%')
PARADIM_UI['s4']['button'].w = Box([PARADIM_UI['s4']['button']])

PARADIM_UI['s4']['status'] = ui.String(description='job status',name='job status',value='')
PARADIM_UI['s4']['status'].disabled = True

PARADIM_UI['s4']['l1'] = VBox([
                               PARADIM_UI['s4']['flfrq'].dd,
                               PARADIM_UI['s4']['qpoints'].dd,
                               PARADIM_UI['s4']['button'].w,
                               PARADIM_UI['s4']['status'].dd
                              
                              ])
PARADIM_UI['s4']['l2'] = VBox([
                               PARADIM_UI['s4']['input'].dd,
                               PARADIM_UI['s4']['output'].dd,
                              ])
PARADIM_UI['s4']['bs'] = HBox([PARADIM_UI['s4']['l1'],PARADIM_UI['s4']['l2']])
PARADIM_UI['s4']['bs'].w = Box([PARADIM_UI['s4']['bs']])
PARADIM_UI['s4']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s4']['bs'].layout = Layout(width='100%', border='1px')

s4_tab0 = ui.Form([PARADIM_UI['s4']['bs']], name = 'Crystal Inputs')
s4_tab1 = ui.Form([PARADIM_UI['s4']['stdin']], name = 'stdin')
s4_tab2 = ui.Form([PARADIM_UI['s4']['stdout']], name = 'stdout')
s4_tab3 = ui.Form([PARADIM_UI['s4']['stderr']], name = 'stderr')


def UpdateStep4( event ):
    global PARADIM_UI, PARADIM
    change = {  
                'fildyn' : 'dyn',
                'flfrc' : PARADIM_UI['s3']['flfrc'].value,        
                'flfrq' : PARADIM_UI['s4']['flfrq'].value,
                'qpoints' : PARADIM_UI['s4']['qpoints'].value,
             }
  
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s4']['input'].value = PARADIM['MATDYN_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s4']['flfrq'].dd.observe(UpdateStep4)
PARADIM_UI['s3']['flfrc'].dd.observe(UpdateStep4)

PARADIM_UI['s4']['display'] = ui.Tab([s4_tab0, s4_tab1, s4_tab2, s4_tab3])
UpdateStep4({'type':'change','name':'value','new':'new'})



##################################################
# Fifth
##################################################


PARADIM_UI['s5'] = {}
PARADIM_UI['s5']['button'] = Button(description='Plot Vibrational Frequencies')
PARADIM_UI['s5']['button'].layout = Layout(width='99%')
PARADIM_UI['s5']['button'].w = Box([PARADIM_UI['s5']['button']])

PARADIM_UI['s5']['input'] = ui.Text( description="inputs", name="inputs", value='''''')
PARADIM_UI['s5']['input'].dd.layout = Layout(width='100%', height='470px')
PARADIM_UI['s5']['input'].disabled = True

PARADIM_UI['s5']['output'] = Output()
PARADIM_UI['s5']['output'].layout = Layout(width='99%')

PARADIM_UI['s5']['l1'] = VBox([
                               PARADIM_UI['s5']['input'].dd,
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

PARADIM_UI['s5']['display'] = ui.Form([PARADIM_UI['s5']['bs']], name = 'Vibrational Frequencies Visualization')