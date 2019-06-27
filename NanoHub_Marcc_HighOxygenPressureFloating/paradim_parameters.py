from ipywidgets import HBox, VBox, Box,  Button, Layout, FloatProgress, Output, Text
import hublib.ui as ui
import random, string
from string import Template
from ipywidgets import Password
from  hublib.ui.formvalue import FormValue
class MyPassword(FormValue):
    def __init__(self, name, value, **kwargs):
        self.dd = Password(value=value)
        FormValue.__init__(self, name, **kwargs)
        
        
PARADIM = {}
PARADIM['SERVER'] = 'gateway2.marcc.jhu.edu'    
PARADIM['TUTORIAL_NAME'] = 'SUMMER_SCHOOL'
PARADIM['USER'] = ''
PARADIM['CODE'] = ''
PARADIM['PWD'] = ''
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
 ibrav = 1,
 a     = 1.00000e+01,
 nat   = 1,
 ntyp  = 1,
 ecutwfc = 30.0,
/
&electrons
 conv_thr = 1.0d-6
/
ATOMIC_SPECIES
${prefix} ${atomic_species} ${upf}
ATOMIC_POSITIONS {crystal}
${prefix}   0.500000   0.500000   0.500000
K_POINTS automatic
2 2 2 1 1 1'''
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

PARADIM['PH_TEMPLATE'] = Template('''Title
&inputph
 prefix = '${prefix}'
 ldisp = .true.
 amass(1) = ${amass},
 fildyn = '${fildyn}',
 nq1 = ${nq1},
 nq2 = ${nq2},
 nq3 = ${nq3},
 tr2_ph = 1.0d-14,
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
PARADIM_UI['s0']['pwd'] = MyPassword(description='MARCC Password', name='MARCC Password',value=PARADIM['PWD'])
PARADIM_UI['s0']['code'] = ui.String(description='GOOGLE AUTH Code', name='GOOGLE AUTH Code',value=PARADIM['CODE'])
PARADIM_UI['s0']['user'] = ui.String(description='Marcc User name', name='Marcc User name',value=PARADIM['USER'])
PARADIM_UI['s0']['folder'] = ui.String(description='Working folder', name='Working folder',value=PARADIM['TUTORIAL_NAME'])

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
PARADIM_UI['s1']['base_mat'] = {
    'O':{ "atomic_species":"15.99940", "upf":"O.pz-n-rrkjus_psl.0.1.UPF" },
    'Ni':{ "atomic_species":"58.69340", "upf":"Ni.pz-n-rrkjus_psl.0.1.UPF" },
    'La':{ "atomic_species":"138.90550", "upf":"La.pz-spfn-rrkjus_psl.1.0.0.UPF" },
}

PARADIM_UI['s1']['commands'] = ui.Text( description="commands",  name="commands", value='''''')
PARADIM_UI['s1']['commands'].dd.layout = Layout(width='99%', height='150px')
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

PARADIM_UI['s1']['button'] = Button(description='Calculate DFT Energy')
PARADIM_UI['s1']['button'].w = Box([PARADIM_UI['s1']['button']])

PARADIM_UI['s1']['job_id'] = Text(placeholder='job ID',value='')
PARADIM_UI['s1']['job_id'].disabled = True
PARADIM_UI['s1']['status'] = Text(placeholder='job status',value='')
PARADIM_UI['s1']['status'].disabled = True
PARADIM_UI['s1']['bc'] = HBox([PARADIM_UI['s1']['button'], PARADIM_UI['s1']['job_id'], PARADIM_UI['s1']['status']])
for mat, param in PARADIM_UI['s1']['base_mat'].items():
    PARADIM_UI['s1']['input_'+mat] = ui.Text( description="O.in",  name="O.in", value='''''')
    PARADIM_UI['s1']['input_'+mat].dd.layout = Layout(width='99%', height='150px')
    PARADIM_UI['s1']['input_'+mat].disabled = True
    PARADIM_UI['s1']['l2_'+mat] = ui.Form([PARADIM_UI['s1']['input_'+mat].dd], name = mat+'.in')

PARADIM_UI['s1']['l2'] = ui.Tab([PARADIM_UI['s1']['l2_'+mat] for mat, params in PARADIM_UI['s1']['base_mat'].items()])

PARADIM_UI['s1']['bs'] = VBox([
    PARADIM_UI['s1']['l2'],
    PARADIM_UI['s1']['bc'],
    PARADIM_UI['s1']['commands'].dd,
])
PARADIM_UI['s1']['bs'].w = Box([PARADIM_UI['s1']['bs']])
PARADIM_UI['s1']['l2'].layout = Layout(width='99%', border='1px')
PARADIM_UI['s1']['bs'].layout = Layout(width='100%', border='1px')

s1_tab0 = ui.Form([PARADIM_UI['s1']['bs']], name = 'Crystal Inputs')
s1_tab1 = ui.Form([PARADIM_UI['s1']['stdin']], name = 'stdin')
s1_tab2 = ui.Form([PARADIM_UI['s1']['stdout']], name = 'stdout')
s1_tab3 = ui.Form([PARADIM_UI['s1']['stderr']], name = 'stderr')


def UpdateStep1( event ):
    global PARADIM_UI, PARADIM
    for mat, params in PARADIM_UI['s1']['base_mat'].items():
        change = {  'k_points':'automatic \n 4 4 4 1 1 1\n', 
                    'prefix':mat,
                    'calculation':'scf',
                    'atomic_species' : params["atomic_species"],
                    'upf': params["upf"],
                 }

        if event['type'] == 'change' and event['name'] == 'value' and event['new']:
            PARADIM_UI['s1']['input_'+mat].value = PARADIM['PW_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s1']['display'] = ui.Tab([s1_tab0, s1_tab1, s1_tab2, s1_tab3])
UpdateStep1({'type':'change','name':'value','new':'new'})

##################################################
# Second
##################################################

PARADIM_UI['s2'] = {}
PARADIM_UI['s2']['amass'] = ui.String(name='amass',description='amass',value='12.0107')
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
                               PARADIM_UI['s2']['amass'].dd,
                               PARADIM_UI['s2']['nq'].dd,
                               PARADIM_UI['s2']['button'].w,
                               PARADIM_UI['s2']['job_id'].dd,
                               PARADIM_UI['s2']['status'].dd,                               
                              ])
PARADIM_UI['s2']['l2'] = VBox([
                               PARADIM_UI['s2']['input'].dd,
                               PARADIM_UI['s2']['commands'].dd,                               
                              ])
PARADIM_UI['s2']['bs'] = HBox([PARADIM_UI['s2']['l1'],PARADIM_UI['s2']['l2']])
PARADIM_UI['s2']['bs'].w = Box([PARADIM_UI['s2']['bs']])
PARADIM_UI['s2']['l2'].layout = Layout(width='100%', border='1px')
PARADIM_UI['s2']['bs'].layout = Layout(width='100%', border='1px')

s2_tab0 = ui.Form([PARADIM_UI['s2']['bs']], name = 'Crystal Inputs')
s2_tab1 = ui.Form([PARADIM_UI['s2']['stdin']], name = 'stdin')
s2_tab2 = ui.Form([PARADIM_UI['s2']['stdout']], name = 'stdout')
s2_tab3 = ui.Form([PARADIM_UI['s2']['stderr']], name = 'stderr')


def UpdateStep2( event ):
    global PARADIM_UI, PARADIM
    change = {  
                'prefix':'diamond',
                'amass' : PARADIM_UI['s2']['amass'].value, 
                'fildyn': 'dyn',
                'nq1': PARADIM_UI['s2']['nq'].value,
                'nq2': PARADIM_UI['s2']['nq'].value,
                'nq3': PARADIM_UI['s2']['nq'].value
    }    
    if event['type'] == 'change' and event['name'] == 'value' and event['new']:
        PARADIM_UI['s2']['input'].value = PARADIM['PH_TEMPLATE'].safe_substitute(change)
    
PARADIM_UI['s2']['amass'].dd.observe(UpdateStep2)
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
    value='''21
 0.500 0.500 0.500
 0.450 0.450 0.450
 0.400 0.400 0.400
 0.350 0.350 0.350
 0.300 0.300 0.300
 0.250 0.250 0.250
 0.200 0.200 0.200
 0.150 0.150 0.150
 0.100 0.100 0.100
 0.050 0.050 0.050
 0.000 0.000 0.000
 0.100 0.000 0.000
 0.200 0.000 0.000
 0.300 0.000 0.000
 0.400 0.000 0.000
 0.500 0.000 0.000
 0.600 0.000 0.000
 0.700 0.000 0.000
 0.800 0.000 0.000
 0.900 0.000 0.000
 1.000 0.000 0.000
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