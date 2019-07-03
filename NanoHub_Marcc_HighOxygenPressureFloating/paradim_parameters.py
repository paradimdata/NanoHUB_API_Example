from ipywidgets import HBox, VBox, Box,  Button, Layout, FloatProgress, Output, Text, Output
import hublib.ui as ui
import random, string
from string import Template
from ipywidgets import Password
from  hublib.ui.formvalue import FormValue
class MyPassword(FormValue):
    def __init__(self, name, value, **kwargs):
        self.dd = Password(value=value)
        FormValue.__init__(self, name, **kwargs)

class Myaction:
    def __init__(self):
        self.callback = None
        
    def on_click(self, on_click):
        self.callback = on_click
    
    def click(self, e):
        if (self.callback != None):
            self.callback(e)
            
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

PARADIM_UI['s2']['button'] = Button(description='Calculate Formation Energy')
PARADIM_UI['s2']['button'].layout = Layout(width='99%')
PARADIM_UI['s2']['button'].w = Box([PARADIM_UI['s2']['button']])
PARADIM_UI['s2']['action'] = Myaction()
PARADIM_UI['s2']['LaNiO3'] = ui.Number(description="LaNiO3", name="LaNiO3", value = -1.35, units="Ha", step=0.1)
PARADIM_UI['s2']['La2O3'] = ui.Number(description="La2O3", name="La2O3", value = -1.485, units="Ha", step=0.1)
PARADIM_UI['s2']['NiO'] = ui.Number(description="NiO", name="NiO", value = -0.475, units="Ha", step=0.1)
PARADIM_UI['s2']['O2'] = ui.Number(description="O2", name="O2", value = -0.23, units="Ha", step=0.1)
PARADIM_UI['s2']['La3Ni2O7'] = ui.Number(description="La3Ni2O7", name="La3Ni2O7", value = -4, units="Ha", step=0.1)
PARADIM_UI['s2']['La4Ni3O10'] = ui.Number(description="La4Ni3O10", name="La4Ni3O10", value = -2.65, units="Ha", step=0.1)
PARADIM_UI['s2']['La2NiO4'] = ui.Number(description="La2NiO4", name="La2NiO4", value = -2, units="Ha", step=0.1)
PARADIM_UI['s2']['energy'] = ui.Number(description="Binding energy", name="Binding energy", value = 0, disabled=True, units="Ha")
                 
    
                 
                 
PARADIM_UI['s2']['button'].on_click(lambda e: PARADIM_UI['s2']['action'].click( e ) )

PARADIM_UI['s2']['l1'] = VBox([
                                PARADIM_UI['s2']['LaNiO3'],
                                PARADIM_UI['s2']['La2O3'],
                                PARADIM_UI['s2']['NiO'],
                                PARADIM_UI['s2']['O2'],
                                PARADIM_UI['s2']['La3Ni2O7'],
                                PARADIM_UI['s2']['La4Ni3O10'],
                                PARADIM_UI['s2']['La2NiO4'],
                                PARADIM_UI['s2']['button'].w,
                                PARADIM_UI['s2']['energy'],
                              ])
PARADIM_UI['s2']['l1'].layout = Layout(width='100%', border='1px')

PARADIM_UI['s2']['display'] = PARADIM_UI['s2']['l1']


##################################################
# Third
##################################################

PARADIM_UI['s3'] = {}

PARADIM_UI['s3']['button'] = Button(description='Calculate Thermodynamics')
PARADIM_UI['s3']['button'].layout = Layout(width='99%')
PARADIM_UI['s3']['button'].w = Box([PARADIM_UI['s3']['button']])
PARADIM_UI['s3']['action'] = Myaction()
PARADIM_UI['s3']['button'].on_click(lambda e: PARADIM_UI['s3']['action'].click( e ) )
PARADIM_UI['s3']['output'] = Output()

PARADIM_UI['s3']['l1'] = VBox([
                               PARADIM_UI['s3']['button'].w,
                              ])
PARADIM_UI['s3']['l1'].layout = Layout(width='100%', border='1px')

PARADIM_UI['s3']['display'] = PARADIM_UI['s3']['l1']



##################################################
# Fourth
##################################################

PARADIM_UI['s4'] = {}

PARADIM_UI['s4']['button'] = Button(description='Calculate Oxygen Pressure')
PARADIM_UI['s4']['button'].layout = Layout(width='99%')
PARADIM_UI['s4']['button'].w = Box([PARADIM_UI['s4']['button']])
PARADIM_UI['s4']['action'] = Myaction()
PARADIM_UI['s4']['button'].on_click(lambda e: PARADIM_UI['s4']['action'].click( e ) )
PARADIM_UI['s4']['output'] = Output()

PARADIM_UI['s4']['l1'] = VBox([
                               PARADIM_UI['s4']['button'].w,
                              ])
PARADIM_UI['s4']['l1'].layout = Layout(width='100%', border='1px')

PARADIM_UI['s4']['display'] = PARADIM_UI['s4']['l1']



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