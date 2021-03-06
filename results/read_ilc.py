import numpy as np
fname = 's4_cmb_ilc.npy'
res_dic = np.load(fname, allow_pickle = 1, encoding = 'latin1').item()
cl_residual = res_dic['cl_residual']
print('\nILC residuals: %s' %(cl_residual.keys()))
for which_spec in cl_residual:
    ilc_nl = cl_residual[which_spec] # ILC residuals
    print('\tSpectra = %s; length = %s' %(which_spec,len(ilc_nl) ))


'''
fg_res_dic = res_dic['fg_res_dic']
print('\nNow residual foreground signals: %s' %(fg_res_dic.keys()))
for which_spec in fg_res_dic:
    print('\tSpectra = %s' %which_spec)
    for gal_res_name in fg_res_dic[which_spec]:
        gal_res_cl = fg_res_dic[which_spec][gal_res_name]
        print('\t\t Signal = %s; length = %s' %(gal_res_name, len(gal_res_cl)))
'''