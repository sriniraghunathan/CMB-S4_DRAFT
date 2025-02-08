import numpy as np, glob, os, sys

batch = True
pgmname = 'get_fisher_withdelensing.py'
ilc_folder = 'ilc_outputs/20210506_with202102designtoolinputforpySM3sims_sedscalingfordust/s4_so_joint_configs/'

#expname_arr = ['s4wide_scaled_sobaseline', 's4wide_scaled_aso', 's4wide_single_chlat_plus_aso', 's4wide']
#expname_arr = ['s4wide_scaled_aso_plus_fulls4scaledsobaseline', 's4wide_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_aso_plus_fulls4scaledsobaseline']
#expname_arr = ['s4wide_single_chlat_plus_2028aso_plus_fulls4scaledsobaseline', 's4wide_single_chlat_plus_2028aso']
#expname_arr = ['s4wide_single_chlat_plus_2028aso']
#expname_arr = ['s4wide_single_chlat_plus_2028aso_plus_fulls4scaledsobaseline']
expname_arr = ['s4deepv3r025']

if batch:
    template = open('batch_jobs/template_slurm.sh')
    #yearval = fname.split('_')[-1].replace('.npy', '').replace('for', '')
    batch_fname = 'batch_jobs/fisher_forecasting.sh'
    opf = open(batch_fname, 'w')
    for line in template:
        opf.writelines('%s\n' %(line.strip()))
for expname in expname_arr:
    searchstr = '%s/%s/*.npy' %(ilc_folder, expname)
    flist = sorted( glob.glob( searchstr ) )
    print('\nExperiment = %s; Total files = %s' %(expname, len( flist ) ))
    for fname in flist:
        opfname = fname.replace('ilc_outputs/', 'results/')
        cmd = 'python3 %s -ilc_fname %s -opfname %s' %(pgmname, fname, opfname)
        if batch:
            opf.writelines('%s\n' %(cmd))
        if not batch:
            os.system(cmd)
            print('\n%s\n' %(cmd))
if batch:
    opf.close()
    cmd = 'sbatch %s' %(batch_fname)
    os.system(cmd)
    print('\n%s\n' %(cmd))
sys.exit()
