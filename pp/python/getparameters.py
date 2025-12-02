# Open copy of input parameterfile
uparamfname = 'used'

#Create copy of parameter file
os.system('cat '+sys.argv[1]+' > '+uparamfname)

# Read from parameterfile into dictionary
dict=pykDict.pykDict()
dict.read_from_file(uparamfname)

# Here we define the variables from the dictionary entries
hpkvd_params      = dict['hpkvd_params']
compile_hpkvd     = dict['compile_hpkvd']
create_filterbank = dict['create_filterbank']
merge_params      = dict['merge_params']
compile_merge     = dict['compile_merge']
map_params        = dict['map_params']
compile_maps      = dict['compile_maps']

batch             = dict['batch']
submit            = dict['submit']

machine           = dict['machine']
submit_command    = dict['submit_command']

seed              = dict['seed']
run_name          = dict['run_name']
short_name        = dict['short_name']
runtype           = dict['runtype']

boxsize           = dict['boxsize']
nmesh             = dict['nmesh']
nbuff             = dict['nbuff']
ntile             = dict['ntile']

tlimit            = dict['tlimit']
largerun          = dict['largerun']
nnodes            = dict['nnodes']
tpnode            = dict['tpnode']
ntasks            = dict['ntasks']
ncpus             = dict['ncpus']
nompth            = dict['nompth']

ievol             = dict['ievol']			
num_redshifts     = dict['num_redshifts']
maximum_redshift  = dict['maximum_redshift']
global_redshift   = dict['global_redshift']

ilpt              = dict['ilpt']

ioutfield         = dict['ioutfield']
ireadfield        = dict['ireadfield']
iwant_field_part  = dict['iwant_field_part']
fielddir          = dict['fielddir']
densfilein        = dict['densfilein']
densfileout       = dict['densfileout']

NonGauss          = dict['NonGauss']
fNL               = dict['fNL']
A_nG              = dict['A_nG']
B_nG              = dict['B_nG']
R_nG              = dict['R_nG']
m_phi             = dict['m_phi']
m_chi             = dict['m_chi']
phi_w             = dict['phi_w']
phi_p             = dict['phi_p']
vev               = dict['vev']
m_tach            = dict['m_tach']
a_e               = dict['a_e']

ntilemerge        = dict['ntilemerge']
ntasksmerge       = dict['ntasksmerge']
iwrap             = dict['iwrap']

Omx               = dict['Omx']
OmB               = dict['OmB']
Omvac             = dict['Omvac']
h                 = dict['h']
ns                = dict['ns']
sigma8            = dict['sigma8']

ioutshear         = dict['ioutshear']	
wsmooth           = dict['wsmooth']
rmax2rs           = dict['rmax2rs']
Rsmooth_max       = dict['Rsmooth_max']

rapi              = dict['rapi']

iforce_strat      = dict['iforce_strat']
ivir_strat        = dict['ivir_strat']
fcoll_3           = dict['fcoll_3']
fcoll_2           = dict['fcoll_2']
fcoll_1           = dict['fcoll_1']
dcrit             = dict['dcrit']

next              = dict['next']
dcore_box         = dict['dcore_box']		
dL_box            = dict['dL_box']	
mlatt             = dict['mlatt']	

cenx              = dict['cenx']
ceny              = dict['ceny']
cenz              = dict['cenz']	
nbuff             = dict['nbuff']


nlx               = dict['nlx']
nly               = dict['nly']
nlz               = dict['nlz']	
n1                = dict['n1']
n2                = dict['n2']
n3                = dict['n3']

filterfile        = dict['filterfile']	

#Merge Parameters
iZeld             = dict['iZeld']
iLexc             = dict['iLexc']
iLmrg             = dict['iLmrg']
iFexc             = dict['iLexc']
iFmrg             = dict['iFmrg']

#Homogeneous ellipsiod collapse table params
TabInterpFile     = dict['TabInterpFile']
TabInterpNx       = dict['TabInterpNx']
TabInterpNy       = dict['TabInterpNy']
TabInterpNz       = dict['TabInterpNz']
TabInterpX1       = dict['TabInterpX1']
TabInterpX2       = dict['TabInterpX2']
TabInterpY1       = dict['TabInterpY1']
TabInterpY2       = dict['TabInterpY2']
TabInterpZ1       = dict['TabInterpZ1']
TabInterpZ2       = dict['TabInterpZ2']

#Mapmaking parameters
maps              = dict['maps']
np_map            = dict['np_map']
nompth_map        = dict['nompth_map']
nnodes_map        = dict['nnodes_map']
tlimit_map        = dict['tlimit_map']
ntasks_map        = dict['ntasks_map']
ppn_map           = dict['ppn_map']
nside_map         = dict['nside_map']
npix_map          = dict['npix_map']  
fov_map           = dict['fov_map']   

zmin_map          = dict['zmin_map']  
zmax_map          = dict['zmax_map'] 

tabfile_map       = dict['tabfile_map']
tabfile_sfr       = dict['tabfile_sfr']
model_map         = dict['model_map']  

scramble_map      = dict['scramble_map']
center_map        = dict['center_map']  
chihview_map      = dict['chihview_map'] 
PSZcut_map        = dict['PSZcut_map']  

ellmax            = dict['ellmax']  

pkfile            = dict['pkfile']

os.system('mv '+uparamfname+' param/'+run_name+'.params')
