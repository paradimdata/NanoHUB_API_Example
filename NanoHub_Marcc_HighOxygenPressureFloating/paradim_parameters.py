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
  calculation = '${calculation}',
  prefix = '${prefix}',
  pseudo_dir = './',
  outdir = './',
  nstep = 300, 
/
&system
  celldm(1)= ${celldm},
  ibrav = 0,
  nat   = ${nat},
  ntyp  = ${ntyp},
  ecutwfc = 30.0,
  occupations = 'smearing',
  smearing    = 'gauss',
  degauss     = 1.00000e-02,
/
&electrons
 conv_thr = 1.0d-6
/
ATOMIC_SPECIES
${atomic_species}
CELL_PARAMETERS ${cell_parameters}
ATOMIC_POSITIONS (crystal)
${atomic_positions}
K_POINTS automatic
4 4 4 1 1 1'''
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


##################################################
# First
##################################################
PARADIM_UI['s1'] = {}
PARADIM_UI['s1']['base_mat'] = {
    'La3Ni2O7' : { 
        "calculation": "scf", 
        'prefix' : 'La3Ni2O7',
        'celldm' : '19.64219137',
        'nat' : '24',
        'ntyp' : '3',
        'units': 2,
        'atomic_species' : '''
La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat= 19.64219137)
    0.990571109  -0.007352894  -0.001439118
    -0.855295511   0.500011253   0.001072431
    0.000700978   0.000202413   0.511712174
''',
        'atomic_positions' : '''
La       0.571601577   0.927918129   0.250195173
La       0.428431896   0.072116175   0.749908240
La       0.928084385   0.571653934   0.250618599
La       0.071942917   0.428372356   0.750117885
La       0.750222120   0.750256434   0.249971124
La       0.249773663   0.249744990   0.750137861
Ni       0.346697828   0.153177209   0.249957563
Ni       0.653312828   0.846846240   0.749994844
Ni       0.153275184   0.346871872   0.249882179
Ni       0.846739908   0.653147589   0.749998378
O        0.095704177   0.904330660   0.999896299
O        0.904419894   0.095382777   0.499780075
O        0.904677229   0.095802784  -0.000229402
O        0.095196218   0.904506635   0.499918698
O        0.595079265   0.404435655   0.500145204
O        0.404785109   0.595791340   0.000107887
O        0.404308189   0.595163649   0.500296535
O        0.595896951   0.404667530   1.000314818
O        0.453554387   0.047209001   0.250885823
O        0.546455316   0.952804650   0.750119600
O        0.750618493   0.750548224   0.750397911
O        0.249393237   0.249495031   0.248719956
O        0.952254964   0.546015429   0.748892343
O        0.047740266   0.453957707   0.250002407        
''',
    },
    'La4Ni3O10' : { 
        "calculation": "scf", 
        'prefix' : 'La4Ni3O10',
        'celldm' : '26.23922530',
        'nat' : '34',
        'ntyp' : '3',
        'units': 2,        
        'atomic_species' : '''
La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat= 26.23922530)
    1.000392219  -0.001102287   0.001081232
    -0.926768440   0.377276475  -0.000312479
    -0.000321458  -0.000683116   0.384407642
''',
        'atomic_positions' : '''
La       0.934136536   0.067540605   0.501022930
La       0.065870535   0.932458933   0.498992948
La       0.565900291   0.432423000   1.000973936
La       0.434108953   0.567580299  -0.000971721
La       0.802409237   0.198070975   0.496351748
La       0.197594564   0.801928298   0.503635487
La       0.697583296   0.301892557   0.995618467
La       0.302419393   0.698104691   0.004415458
Ni       0.500003242   0.500002882   0.500001654
Ni       0.000000322   0.000003150   1.000009681
Ni       0.639393238   0.359865012   0.498833997
Ni       0.360613151   0.640131807   0.501176244
Ni       0.860652583   0.140139606   0.999607725
Ni       0.139354357   0.859857718   0.000400586
O        0.749903852   0.743720384   0.753616610
O        0.750174638   0.756382218   0.253491027
O        0.250102522   0.256280441   0.246395898
O        0.249832257   0.243620092   0.746518971
O        0.570271976   0.429269280   0.477309534
O        0.429737843   0.570739366   0.522764904
O        0.929845906   0.070807203   0.976966090
O        0.070157809   0.929194405   0.022982533
O        0.886487401   0.613200066   0.749418070
O        0.613838949   0.887062258   0.250361079
O        0.113516366   0.386782446   0.250588056
O        0.386170005   0.112927221   0.749643345
O        0.714857854   0.284587622   0.514071763
O        0.285141988   0.715403617   0.485893950
O        0.785259431   0.215477474   0.014814318
O        0.214746803   0.784519670   0.985241991
O        0.890680275   0.608056658   0.249266662
O        0.609124855   0.891784866   0.750420081
O        0.109317509   0.391946507   0.750744049
O        0.390882060   0.108228674   0.249591931       
''',
    }, 
    'LaNiO3' : { 
        "calculation": "scf", 
        'prefix' : 'LaNiO3',
        'celldm' : '7.28886267',
        'nat' : '5',
        'units': 1,        
        'ntyp' : '3',
        'atomic_species' : '''
La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat=  7.28886267)
    0.973022005   0.000000000   0.000000000
    0.000000000   0.973022005   0.000000000
    0.000000000   0.000000000   0.973022005
''',
        'atomic_positions' : '''
La       0.500000000   0.500000000   0.500000000
Ni       0.000000000   0.000000000   0.000000000
O        0.000000000   0.000000000   0.500000000
O        0.500000000   0.000000000   0.000000000
O        0.000000000   0.500000000   0.000000000  
''',
    },
    'La2O3' : { 
        "calculation": "scf", 
        'prefix' : 'La2O3',
        'celldm' : '26.23922530',
        'nat' : '5',
        'units': 1,        
        'ntyp' : '2',
        'atomic_species' : '''
La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat=  7.19040794)
    0.987661282   0.000000618  -0.000001548
    -0.493830106   0.855340070   0.000001548
    -0.000004065   0.000002347   2.730934463
''',
        'atomic_positions' : '''
La       0.333333076   0.666666924   0.647244330
La       0.666666924   0.333333076   0.352755670
O        0.333332764   0.666667236   0.286911985
O        0.666667236   0.333332764   0.713088015
O        0.000000000   0.000000000   0.500000000      
''',
    },
    'NiO' : { 
        "calculation": "scf", 
        'prefix' : 'NiO',
        'celldm' : '8.19194389',
        'nat' : '2',
        'units': 1,        
        'ntyp' : '2',
        'atomic_species' : '''
Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat=  8.19194389)
    -0.466966682  -0.000000000   0.466966682
    0.000000000   0.466966682   0.466966682
    -0.466966682   0.466966682  -0.000000000
''',
        'atomic_positions' : '''
Ni       0.000000000  -0.000000000  -0.000000000
O        0.500000000   0.500000000   0.500000000     
''',
    },
    'O2' : { 
        "calculation": "scf", 
        'prefix' : 'O2',
        'celldm' : '23.68015817',
        'nat' : '2',
        'units': 1,        
        'ntyp' : '1',
        'atomic_species' : '''
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat= 23.68015817)
    1.000088960   0.000000000   0.000000000
    0.000000000   1.000088960   0.000000000
    0.000000000   0.000000000   0.769820460
''',
        'atomic_positions' : '''
O        0.500000000   0.500000000   0.937318663
O        0.500000000   0.500000000   0.062681337     
''',
    },    
    
    'La2NiO4' : { 
        "calculation": "scf", 
        'prefix' : 'La2NiO4',
        'celldm' : '13.41856732',
        'nat' : '7',
        'units': 1,        
        'ntyp' : '3',
        'atomic_species' : '''
La    138.90550  La.pz-spfn-rrkjus_psl.1.0.0.UPF
Ni     58.69340  Ni.pz-n-rrkjus_psl.0.1.UPF
O      15.99940  O.pz-n-rrkjus_psl.0.1.UPF
''',
        'cell_parameters' : '''(alat= 13.41856732)
    0.981971436  -0.003972772   0.000001791
    0.709373515   0.679023413   0.000001933
    -0.845673282  -0.337525566   0.367691005
''',
        'atomic_positions' : '''
La       0.637550915   0.637550903   0.000000101
La       0.362449085   0.362449097  -0.000000101
O        0.818296061   0.818295711   0.000000189
O        0.181703939   0.181704289  -0.000000189
Ni       0.000000000   0.000000000   0.000000000
O        0.000002000   0.500002000   0.500003000
O        0.500002000   0.000002000   0.500003000
''',
    },
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
PARADIM_UI['s1']['action'] = Myaction()
PARADIM_UI['s1']['button'].on_click(lambda e: PARADIM_UI['s1']['action'].click( e ) )


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
        change = { 
                    'prefix':params["prefix"],
                    'calculation':params["calculation"],
                    'celldm':params['celldm'],
                    'nat':params['nat'],
                    'ntyp':params['ntyp'],
                    'atomic_species':params['atomic_species'],
                    'cell_parameters':params['cell_parameters'],   
                    'atomic_positions':params['atomic_positions'],
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
PARADIM_UI['s2']['LaNiO3'] = ui.Number(description="LaNiO3", name="LaNiO3 ~ -325.81", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['La2O3'] = ui.Number(description="La2O3", name="La2O3 ~ -350.69", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['NiO'] = ui.Number(description="NiO", name="NiO ~-133.62", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['O2'] = ui.Number(description="O2", name="O2 ~ -66.67", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['La3Ni2O7'] = ui.Number(description="La3Ni2O7", name="La3Ni2O7 ~ -810.27", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['La4Ni3O10'] = ui.Number(description="La4Ni3O10", name="La4Ni3O10 ~ -1136.11", value = 0, units="Ry", step=0.1)
PARADIM_UI['s2']['La2NiO4'] = ui.Number(description="La2NiO4", name="La2NiO4 ~ -484.42", value = 0, units="Ry", step=0.1)

PARADIM_UI['s2']['DE_LaNiO3'] = ui.Number(description="Binding energy LaNiO3", name="Binding energy LaNiO3", value = 0, disabled=True, units="Ha")
PARADIM_UI['s2']['DE_La3Ni2O7'] = ui.Number(description="Binding energy La3Ni2O7", name="Binding energy La3Ni2O7", value = 0, disabled=True, units="Ha")
PARADIM_UI['s2']['DE_La4Ni3O10'] = ui.Number(description="Binding energy La4Ni3O10", name="Binding energy La4Ni3O10", value = 0, disabled=True, units="Ha")
PARADIM_UI['s2']['DE_La2NiO4'] = ui.Number(description="Binding energy La2NiO4", name="Binding energy La2NiO4", value = 0, disabled=True, units="Ha")
                 
    
                 
                 
PARADIM_UI['s2']['button'].on_click(lambda e: PARADIM_UI['s2']['action'].click( e ) )

PARADIM_UI['s2']['l1'] = HBox([
                            VBox([    
                                ui.Form([
                                    PARADIM_UI['s2']['LaNiO3'],
                                    PARADIM_UI['s2']['La2O3'],
                                    PARADIM_UI['s2']['NiO'],
                                    PARADIM_UI['s2']['O2'],
                                    PARADIM_UI['s2']['La3Ni2O7'],
                                    PARADIM_UI['s2']['La4Ni3O10'],
                                    PARADIM_UI['s2']['La2NiO4'],
                                ], name="Total Energies"),
                            ], layout = Layout(width='100%', padding="2px")),
                            VBox([
                                PARADIM_UI['s2']['DE_LaNiO3'],
                                PARADIM_UI['s2']['DE_La3Ni2O7'],
                                PARADIM_UI['s2']['DE_La4Ni3O10'],
                                PARADIM_UI['s2']['DE_La2NiO4'],
                                PARADIM_UI['s2']['button'].w
                            ], layout = Layout(width='100%', padding="2px"))
                        ], layout = Layout(width='100%', border='1px'))

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
