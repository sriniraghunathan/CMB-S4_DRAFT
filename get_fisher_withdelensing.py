import sys, scipy, pickle, numpy, argparse, os, copy, time
#import cambWrapTools
import classWrapTools
import fisherTools
#from mpi4py import MPI

parser = argparse.ArgumentParser(description='')
parser.add_argument('-lmin', dest='lmin', action='store', help='lmin', type=int, default=30)
parser.add_argument('-lmax', dest='lmax', action='store', help='lmax', type=int, default=5000)
parser.add_argument('-lmaxTT', dest='lmaxTT', action='store', help='lmaxTT for TT Fisher', type=int, default=5000)
parser.add_argument('-lmax_T_lensing', dest='lmax_T_lensing', action='store', help='CMB lmax for T-based lensing', type=int, default=3000)
parser.add_argument('-lmax_P_lensing', dest='lmax_P_lensing', action='store', help='CMB lmax for P-based lensing', type=int, default=3000)
parser.add_argument('-lbuffer', dest='lbuffer', action='store', help='lbuffer for non-Gaussian stuffs', type=int, default=0)
parser.add_argument('-delta_l_max', dest='delta_l_max', action='store', help='delta_l_max for non-Gaussian stuffs', type=int, default=0)
parser.add_argument('-ilc_fname', dest='ilc_fname', action='store', help='ilc_fname', type=str, default='ilc_outputs/s4deepv3r025_ilc_galaxy1_27-39-93-145-225-278_TT-EE_galmask5_AZ_for7years.npy')

args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)

start = time.time()
#MPI
'''
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(rank, size)
'''
rank = 0
size = 1

###  Set of experiments  ###
# Results will be indexed by experiment number, starting from 0
#expNames = list(range(3))
expNames = [0] #srini - just 1 experiment
nExps = len(expNames)
'''
lmax = 5000
lmaxTT = 5000
lmin = 30
#noiseLevels = numpy.arange(0.5, 10.5, 0.5)
#beamSizeArcmin = 1.0

lbuffer = 1500
'''
lmax_calc = lmax+lbuffer

#expNamesThisNode = numpy.array_split(numpy.asarray(expNames), size)[rank]

# Directory where CLASS_delens is located
classExecDir = './CLASS_delens/'
# Directory where you would like the output
classDataDir = './CLASS_delens/'
outputDir = classDataDir + 'results/'

classDataDirThisNode = classDataDir + 'data/Node_' + str(rank) + '/'
# Base name to use for all output files
fileBase = 'fisher_8p'
fileBaseThisNode = fileBase + '_' + str(rank)

if not os.path.exists(classDataDirThisNode):
    os.makedirs(classDataDirThisNode)
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

# Spectrum types and polarizations to include
spectrumTypes = ['unlensed', 'lensed', 'delensed', 'lensing']
polCombs = ['cl_TT', 'cl_TE', 'cl_EE', 'cl_dd']

#######################################################################################3
#LOAD PARAMS AND GET POWER SPECTRA

#Fiducial values and step sizes taken from arXiv:1509.07471 Allison et al
cosmoFid = {'omega_c_h2':0.1197, \
                'omega_b_h2': 0.0222, \
                'N_eff': 3.046, \
                'A_s' : 2.196e-9, \
                'n_s' : 0.9655,\
                'tau' : 0.06, \
                #'H0' : 67.5, \
                'theta_s' : 0.010409, \
                #'Yhe' : 0.25, \
                #'r'   : 0.01, \
                'mnu' : 0.06}
#cosmoFid['n_t'] = - cosmoFid['r'] / 8.0 * (2.0 - cosmoFid['n_s'] - cosmoFid['r'] / 8.0)

stepSizes = {'omega_c_h2':0.0030, \
                'omega_b_h2': 0.0008, \
                'N_eff': .080, \
                'A_s' : 0.1e-9, \
                'n_s' : 0.010,\
                'tau' : 0.020, \
                'H0' : 1.2, \
                'theta_s' : 0.000050, \
                'mnu' : 0.02, \
                #'r'   : 0.001, \
                #'n_t' : cosmoFid['n_t'], \
                'Yhe' : 0.0048}

cosmoParams = list(cosmoFid.keys())
#srini
##delta_l_max = 2000 
ell = numpy.arange(2,lmax_calc+1+delta_l_max)

# Mask the \ells you do not want included in lensing reconstruction
# Keys can be added as e.g. 'lmin_T', 'lmax_T', etc.
reconstructionMask = dict()
reconstructionMask['lmax_T'] = lmax_T_lensing
#srini
reconstructionMask['lmax_E'] = lmax_P_lensing
reconstructionMask['lmax_B'] = lmax_P_lensing

extra_params = dict()
#extra_params['delensing_verbose'] = 3
#extra_params['output_spectra_noise'] = 'no'
#extra_params['write warnings'] = 'y'
extra_params['delta_l_max'] = delta_l_max

# Specify \ells to keep when performing Fisher matrix sum
ellsToUse = {'cl_TT': [lmin, lmaxTT], 'cl_TE': [lmin, lmax], 'cl_EE': [lmin, lmax], 'cl_dd': [2, lmax]}
ellsToUseNG = {'cl_TT': [lmin, lmaxTT], 'cl_TE': [lmin, lmax], 'cl_EE': [lmin, lmax], 'cl_dd': [2, lmax], 'lmaxCov': lmax_calc}

cmbNoiseSpectra = dict()
deflectionNoises = dict()
paramDerivs = dict()
powersFid = dict()
invCovDotParamDerivs_delensed = dict()
invCovDotParamDerivs_lensed = dict()
paramDerivStack_delensed = dict()
paramDerivStack_lensed = dict()
fisherGaussian = dict()
fisherNonGaussian_delensed = dict()
fisherNonGaussian_lensed = dict()

# Flags for whether to include NonGaussian covariances, and derivatives wrt unlensed spectra
doNonGaussian = False #True
includeUnlensedSpectraDerivatives = False #True

# Calculations begin

### Assign task of computing lensed NG covariance to last node       ###
### This is chosen because last node sometimes has fewer experiments ###
if doNonGaussian is True:
    if rank == size-1:

        if includeUnlensedSpectraDerivatives:
            dCldCLd_lensed, dCldCLu_lensed = classWrapTools.class_generate_data(cosmoFid,
                                         cmbNoise = None, \
                                         deflectionNoise = None, \
                                         extraParams = extra_params, \
                                         rootName = fileBaseThisNode, \
                                         lmax = lmax_calc, \
                                         calculateDerivatives = 'lensed', \
                                         includeUnlensedSpectraDerivatives = includeUnlensedSpectraDerivatives,
                                         classExecDir = classExecDir,
                                         classDataDir = classDataDirThisNode)
        else:
            dCldCLd_lensed = classWrapTools.class_generate_data(cosmoFid,
                                         cmbNoise = None, \
                                         deflectionNoise = None, \
                                         extraParams = extra_params, \
                                         rootName = fileBaseThisNode, \
                                         lmax = lmax_calc, \
                                         calculateDerivatives = 'lensed', \
                                         classExecDir = classExecDir,
                                         classDataDir = classDataDirThisNode)
            dCldCLu_lensed = None

        print('Successfully computed derivatives')
    else:
        dCldCLd_lensed = None

planckNoise = fisherTools.getPlanckInvVarNoise(ells = ell, includePol = True)
'''
# for importing DRAFT results
#galaxy = ['0', '1', '1']
#galmask = ['', '_galmask2_AZ', '_galmask5_AZ']
include_gal = 1
#galmask = [0, 1, 2, 3, 4, 5]
galmask = 5
totobsyears = 7
'''

'''
for k in expNamesThisNode:
    expName = expNames[k]
'''
#srini
k=0
expName = 0
if (1):
    print('Node ' + str(rank) + ' working on experiment ' + str(expName))

    # DRAFT import begins
    #fname = '/users/ctrendafilova/scratch/DRAFT/results/20200701/s4like_mask_v2/TT-EE/baseline/s4wide_ilc_galaxy1_27-39-93-145-225-278_TT-EE' + str(galmask[k]) + '.npy'
    #ilc_folder = '/Users/sraghunathan/Research/SPTpol/analysis/git/DRAFT/results/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4deepv3r025/splat_minobsel30_galcuts_mask/TT-EE/baseline'
    '''
    ilc_folder = 'ilc_outputs/'

    if not include_gal:
        fname = '%s/s4deepv3r025_ilc_galaxy%s_27-39-93-145-225-278_TT-EE_for%syears.npy' %(ilc_folder, include_gal, totobsyears)
    else:   
        fname = '%s/s4deepv3r025_ilc_galaxy%s_27-39-93-145-225-278_TT-EE_galmask%s_AZ_for%syears.npy' %(ilc_folder, include_gal, galmask, totobsyears)
    '''
    fname = ilc_fname

    DRAFT_ILC_dic = numpy.load(fname, allow_pickle = 1).item()

    noise_ell = DRAFT_ILC_dic['el'][2:]
    tempNoise = DRAFT_ILC_dic['cl_residual']['TT'][2:]
    polNoise = DRAFT_ILC_dic['cl_residual']['EE'][2:]

    #srini
    ##tempNoise[:20] += 1.  ## Lowest ells give zero from DRAFT, Planck will dominate on these scales
    ##polNoise[:20] += 1.
    tempNoise[noise_ell<11] += 1.  ## Lowest ells give zero from DRAFT, Planck will dominate on these scales
    polNoise[noise_ell<11] += 1.

    '''
    #srini
    ##NlTT = numpy.pad(tempNoise, (0,2000), 'linear_ramp', end_values = (1.,)) # Add higher ells for convenience (not used in calculation)
    ##NlEE = numpy.pad(polNoise, (0,2000), 'linear_ramp', end_values = (1.,)) # Add higher ells for convenience (not used in calculation)
    ell_pad = 5000
    NlTT = numpy.pad(tempNoise, (0,ell_pad), 'linear_ramp', end_values = (1.,)) # Add higher ells for convenience (not used in calculation)
    NlEE = numpy.pad(polNoise, (0,ell_pad), 'linear_ramp', end_values = (1.,)) # Add higher ells for convenience (not used in calculation)
    '''
    NlTT = numpy.interp(ell, noise_ell, tempNoise)
    NlEE = numpy.interp(ell, noise_ell, polNoise)

    cmbNoiseSpectra[k] = classWrapTools.noiseSpectra(l = ell, noiseLevelT = 0, beamArcmin = 0)

    cmbNoiseSpectra[k]['cl_TT'] += NlTT[:len(ell)]
    cmbNoiseSpectra[k]['cl_EE'] += NlEE[:len(ell)]
    cmbNoiseSpectra[k]['cl_BB'] += NlEE[:len(ell)]

    cmbNoiseSpectra[k]['dl_TT'] = ell*(ell+1.) / (2.*numpy.pi) * cmbNoiseSpectra[k]['cl_TT']
    cmbNoiseSpectra[k]['dl_EE'] = ell*(ell+1.) / (2.*numpy.pi) * cmbNoiseSpectra[k]['cl_EE']
    cmbNoiseSpectra[k]['dl_BB'] = ell*(ell+1.) / (2.*numpy.pi) * cmbNoiseSpectra[k]['cl_EE']

    #### Combine Planck and S4 with inverse variance weighting
    noisePolCombs = cmbNoiseSpectra[k].keys()
    TotalCombinedNoise = fisherTools.onedDict(noisePolCombs)
    for pc, polComb in enumerate(noisePolCombs):
        runningTotal = numpy.zeros(len(ell))
        runningTotal += 1/cmbNoiseSpectra[k][polComb]
        if k!=2: ###Exclude Planck for Mask2(b)
            runningTotal += 1/planckNoise[polComb]
        TotalCombinedNoise[polComb] = 1/runningTotal
    #the ells got messed up as part of the process ... fix this.
    TotalCombinedNoise['l'] = ell

    cmbNoiseSpectra[k] = TotalCombinedNoise.copy()
    # copy noise curve constructed from DRAFT + Planck and proceed as usual
    powersFid[k], deflectionNoises[k] = classWrapTools.class_generate_data(cosmoFid,
                                         cmbNoise = cmbNoiseSpectra[k],
                                         extraParams = extra_params,
                                         rootName = fileBaseThisNode,
                                         lmax = lmax_calc,
                                         classExecDir = classExecDir,
                                         classDataDir = classDataDirThisNode,
                                         reconstructionMask = reconstructionMask)

    paramDerivs[k] = fisherTools.getPowerDerivWithParams(cosmoFid = cosmoFid, \
                            extraParams = extra_params, \
                            stepSizes = stepSizes, \
                            polCombs = polCombs, \
                            cmbNoiseSpectraK = cmbNoiseSpectra[k], \
                            deflectionNoisesK = deflectionNoises[k], \
                            useClass = True, \
                            lmax = lmax_calc, \
                            fileNameBase = fileBaseThisNode, \
                            classExecDir = classExecDir, \
                            classDataDir = classDataDirThisNode)

    fisherGaussian[k] = fisherTools.getGaussianCMBFisher(powersFid = powersFid[k], \
                            paramDerivs = paramDerivs[k], \
                            cmbNoiseSpectra = cmbNoiseSpectra[k], \
                            deflectionNoises = deflectionNoises[k], \
                            cosmoParams = cosmoParams, \
                            spectrumTypes = ['unlensed', 'lensed', 'delensed'], \
                            polCombsToUse = polCombs, \
                            ellsToUse = ellsToUse)

    if doNonGaussian:

        ### Overwrite dCldCLd_delensed for each experiment to save memory ###

        if includeUnlensedSpectraDerivatives:
            dCldCLd_delensed, dCldCLu_delensed = classWrapTools.class_generate_data(cosmoFid,
                                                 cmbNoise = cmbNoiseSpectra[k], \
                                                 deflectionNoise = deflectionNoises[k], \
                                                 extraParams = extra_params, \
                                                 rootName = fileBaseThisNode, \
                                                 lmax = lmax_calc, \
                                                 calculateDerivatives = 'delensed', \
                                                 includeUnlensedSpectraDerivatives = includeUnlensedSpectraDerivatives,
                                                 classExecDir = classExecDir,
                                                 classDataDir = classDataDirThisNode)
        else:
            dCldCLd_delensed = classWrapTools.class_generate_data(cosmoFid,
                                                 cmbNoise = cmbNoiseSpectra[k], \
                                                 deflectionNoise = deflectionNoises[k], \
                                                 extraParams = extra_params, \
                                                 rootName = fileBaseThisNode, \
                                                 lmax = lmax_calc, \
                                                 calculateDerivatives = 'delensed', \
                                                 classExecDir = classExecDir,
                                                 classDataDir = classDataDirThisNode)
            dCldCLu_delensed = None


        invCovDotParamDerivs_delensed[k], paramDerivStack_delensed[k] = fisherTools.choleskyInvCovDotParamDerivsNG(powersFid = powersFid[k], \
                                    cmbNoiseSpectra = cmbNoiseSpectra[k], \
                                    deflectionNoiseSpectra = deflectionNoises[k], \
                                    dCldCLd = dCldCLd_delensed,
                                    paramDerivs = paramDerivs[k], \
                                    cosmoParams = cosmoParams, \
                                    dCldCLu = dCldCLu_delensed, \
                                    ellsToUse = ellsToUseNG, \
                                    polCombsToUse = polCombs, \
                                    spectrumType = 'delensed')


        if rank != size-1 and dCldCLd_lensed is None:
            classDataDirLastNode = classDataDir + 'data/Node_' + str(size-1) + '/'
            fileBaseLastNode = fileBase + '_' + str(size-1)

            dCldCLd_lensed = classWrapTools.loadLensingDerivatives(rootName = fileBaseLastNode,
                                                                   classDataDir = classDataDirLastNode,
                                                                   dervtype = 'lensed')


            dCldCLu_lensed = None
            if includeUnlensedSpectraDerivatives:
                dCldCLu_lensed = classWrapTools.loadUnlensedSpectraDerivatives(rootName = fileBaseLastNode,
                                                                   classDataDir = classDataDirLastNode,
                                                                   dervtype = 'lensed')


        invCovDotParamDerivs_lensed[k], paramDerivStack_lensed[k] = fisherTools.choleskyInvCovDotParamDerivsNG(powersFid = powersFid[k], \
                                    cmbNoiseSpectra = cmbNoiseSpectra[k], \
                                    deflectionNoiseSpectra = deflectionNoises[k], \
                                    dCldCLd = dCldCLd_lensed,
                                    paramDerivs = paramDerivs[k], \
                                    cosmoParams = cosmoParams, \
                                    dCldCLu = dCldCLu_lensed,
                                    ellsToUse = ellsToUseNG, \
                                    polCombsToUse = polCombs, \
                                    spectrumType = 'lensed')

        fisherNonGaussian_delensed[k] = fisherTools.getNonGaussianCMBFisher(invCovDotParamDerivs = invCovDotParamDerivs_delensed[k], \
                                    paramDerivStack = paramDerivStack_delensed[k], \
                                    cosmoParams = cosmoParams)

        fisherNonGaussian_lensed[k] = fisherTools.getNonGaussianCMBFisher(invCovDotParamDerivs = invCovDotParamDerivs_lensed[k], \
                                    paramDerivStack = paramDerivStack_lensed[k], \
                                    cosmoParams = cosmoParams)

print('Node ' + str(rank) + ' finished all experiments')

forecastData = {'cmbNoiseSpectra' : cmbNoiseSpectra,
                'powersFid' : powersFid,
                'paramDerivs': paramDerivs,
                'fisherGaussian': fisherGaussian,
                'deflectionNoises' : deflectionNoises}
if doNonGaussian:
    forecastData['invCovDotParamDerivs_delensed'] = invCovDotParamDerivs_delensed
    forecastData['paramDerivStack_delensed'] = paramDerivStack_delensed
    forecastData['invCovDotParamDerivs_lensed'] = invCovDotParamDerivs_lensed
    forecastData['paramDerivStack_lensed'] = paramDerivStack_lensed
    forecastData['fisherNonGaussian_delensed'] = fisherNonGaussian_delensed
    forecastData['fisherNonGaussian_lensed'] = fisherNonGaussian_lensed

print('Node ' + str(rank) + ' saving data')

filename = classDataDirThisNode + fileBaseThisNode + '.pkl'
delensedOutput = open(filename, 'wb')
pickle.dump(forecastData, delensedOutput, -1)
delensedOutput.close()
print('Node ' + str(rank) + ' saving data complete')
end = time.time()
print('\n\n\ttotal time take = %.3f minutes' %((end-start)/60.))
sys.exit()

comm.Barrier()

if rank==0:
    print('Node ' + str(rank) + ' collecting data')
    for irank in range(1,size):
        print('Getting data from node ' + str(irank))
        filename = classDataDir + 'data/Node_' + str(irank) + '/' + fileBase + '_' + str(irank) + '.pkl'
        nodeData = open(filename, 'rb')
        nodeForecastData = pickle.load(nodeData)
        nodeData.close()
        for key in list(forecastData.keys()):
            forecastData[key].update(nodeForecastData[key])

    print('Node ' + str(rank) + ' reading script')
    f = open(os.path.abspath(__file__), 'r')
    script_text = f.read()
    f.close()

    forecastData['script_text'] = script_text

    forecastData['cosmoFid'] = cosmoFid
    forecastData['cosmoParams'] = cosmoParams

    print('Node ' + str(rank) + ' saving collected data')
    filename = outputDir + fileBase + '.pkl'
    delensedOutput = open(filename, 'wb')
    pickle.dump(forecastData, delensedOutput, -1)
    delensedOutput.close()
