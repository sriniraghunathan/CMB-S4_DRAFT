import numpy as np, sys, scipy as sc, os
#from camb import model, initialpower
#from pylab import *
from scipy import linalg
from copy import deepcopy


########################################################################################################################
def get_nldic(nlfile, els):
    dic = np.load(nlfile, allow_pickle = 1, encoding = 'latin1').item()
    el_nl, cl_residual = dic['el'], dic['cl_residual']

    if 'T' in cl_residual:
        nl_TT, nl_EE = cl_residual['T'], cl_residual['P']
        nl_TT = np.interp(els, el_nl, nl_TT)
        nl_EE = np.interp(els, el_nl, nl_EE)
        nl_TE = None
    else:
        nl_TT, nl_EE = cl_residual['TT'], cl_residual['EE']
        if 'TE' in cl_residual:
            nl_TE = cl_residual['TE']
        else:
            nl_TE = None
        el_nl = np.arange(len(nl_TT))
        nl_TT = np.interp(els, el_nl, nl_TT)
        nl_EE = np.interp(els, el_nl, nl_EE)
        if nl_TE is not None:
            nl_TE = np.interp(els, el_nl, nl_TE)
        else:
            nl_TE = np.zeros(len(els))

    nldic = {}
    nldic['TT'] = nl_TT
    nldic['EE'] = nl_EE
    nldic['TE'] = nl_TE

    return nldic, dic['fsky_val']

########################################################################################################################

def fn_ini_param_dic(fpath = 'params/params_planck_r_0.0_2015_cosmo_lensed_LSS.txt'):
    """
    read params file and initialise cosmology
    """
    try:
        params = np.recfromtxt(fpath, delimiter = '=', encoding = 'utf-8')
    except:
        params = np.recfromtxt(fpath, delimiter = '=')
    param_dict = {}
    for rec in params:
        val = rec[1].strip()##.decode("utf-8")
        try:
            if val.find('.')>-1:
                val = float(val)
            else:
                val = int(val)
        except:
            val = str(val)

        if val == 'None':
            val = None
        paramname = rec[0].strip()#.decode("utf-8")
        param_dict[paramname] = val

    return param_dict

########################################################################################################################

def fn_delta_Cl(els, cl_dic, nl_dic, fsky):
    delta_cl_dic = {}
    for XX in cl_dic:
        if XX == 'TT':
            nl = nl_dic['TT']
        elif XX == 'EE' or XX == 'BB':
            nl = nl_dic['EE']
        elif XX == 'TE':
            nl = nl_dic['TE']
            if (1):
                #print('\n\n\n\t\t\t\tnulling galaxy TE\n\n\n')
                nl = np.copy(nl) * 0.

        cl = cl_dic[XX]
        delta_cl_dic[XX] = np.sqrt(2./ (2.*els + 1.) / fsky ) * (cl + nl)

    return delta_cl_dic

########################################################################################################################

def get_cov(TT, EE, TE, PP, TP, EP):

    C = np.zeros( (3,3) ) #TT, EE, PP
    C[0,0] = TT
    C[1,1] = EE
    C[0,1] = C[1,0] = TE

    C[2,2] = PP
    C[0,2] = C[2,0] = TP
    C[1,2] = C[2,1] = 0. ##EP

    return np.mat( C )

########################################################################################################################

def get_fisher_mat(els, cl_deriv_dic, delta_cl_dic, params, pspectra_to_use, min_l_temp = None, max_l_temp = None, min_l_pol = None, max_l_pol = None):

    if min_l_temp is None: min_l_temp = 0
    if max_l_temp is None: max_l_temp = 10000

    if min_l_pol is None: min_l_pol = 0
    if max_l_pol is None: max_l_pol = 10000

    npar = len(params)
    F = np.zeros([npar,npar])
    #els = np.arange( len( delta_cl_dic.values()[0] ) )

    with_lensing = 0
    if 'PP' in pspectra_to_use:
        with_lensing = 1

    all_pspectra_to_use = []
    for tmp in pspectra_to_use:
        if isinstance(tmp, list):      
            all_pspectra_to_use.extend(tmp)
        else:
            all_pspectra_to_use.append(tmp)

    for lcntr, l in enumerate( els ):

        TT, EE, TE = 0., 0., 0.
        Tphi = Ephi = PP = 0.
        if 'TT' in delta_cl_dic:
            TT = delta_cl_dic['TT'][lcntr]
        if 'EE' in delta_cl_dic:
            EE = delta_cl_dic['EE'][lcntr]
        if 'TE' in delta_cl_dic:
            TE = delta_cl_dic['TE'][lcntr]
        if with_lensing:
            Tphi, Ephi, PP = delta_cl_dic['Tphi'][lcntr], delta_cl_dic['Ephi'][lcntr], delta_cl_dic['PP'][lcntr]

        ##############################################################################
        #null unused fields
        null_TT, null_EE, null_TE = 0, 0, 0
        if l<min_l_temp or l>max_l_temp:
            null_TT = 1
        if l<min_l_pol or l>max_l_pol: 
            null_EE = 1
            null_TE = 1
        null_PP = 0 #Lensing noise curves already have pretty large noise outside desired L range

        if 'TT' not in all_pspectra_to_use:
            null_TT = 1
        if 'EE' not in all_pspectra_to_use:
            null_EE = 1
        if 'TE' not in all_pspectra_to_use:
            #if 'TT' not in pspectra_to_use and 'EE' not in pspectra_to_use:
            #    null_TE = 1
            if 'TT' in pspectra_to_use and 'EE' in pspectra_to_use:
                null_TE = 0
            else:
                null_TE = 1
        if ['PP'] not in pspectra_to_use:
            null_PP = 1
        if ['TT', 'EE', 'TE'] in pspectra_to_use:
            null_TT = 0
            null_EE = 0
            null_TE = 0

        #nulling unwanted fields
        if null_TT and null_TE: TT = 0
        if null_EE and null_TE: EE = 0
        #if null_TE and (null_TT and null_EE): TE = 0
        if null_TE: 
            if not null_TT and not null_EE:
                pass
            else:
                TE = 0
        if null_PP: PP = Tphi = EPhi = 0
        if null_TT: Tphi = 0
        if null_EE: Ephi = 0

        if (null_TT and null_EE and null_TE and null_PP): continue
        ##############################################################################
        #get covariance matrix and its inverse
        COV_mat_l = get_cov(TT, EE, TE, PP, Tphi, Ephi)
        inv_COV_mat_l = linalg.pinv2(COV_mat_l)

        ##############################################################################
        #get the parameter combinations
        param_combinations = []
        for pcnt,p in enumerate(params):
            for pcnt2,p2 in enumerate(params):
                ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
                param_combinations.append([p,p2, pcnt, pcnt2])

        ##############################################################################

        for (p,p2, pcnt, pcnt2) in param_combinations:


            TT_der1, EE_der1, TE_der1 = 0., 0., 0.
            TT_der2, EE_der2, TE_der2 = 0., 0., 0.

            if 'TT' in cl_deriv_dic[p]:
                TT_der1 = cl_deriv_dic[p]['TT'][lcntr]
                TT_der2 = cl_deriv_dic[p2]['TT'][lcntr]
            if 'EE' in cl_deriv_dic[p]:
                EE_der1 = cl_deriv_dic[p]['EE'][lcntr]
                EE_der2 = cl_deriv_dic[p2]['EE'][lcntr]
            if 'TE' in cl_deriv_dic[p]:
                TE_der1 = cl_deriv_dic[p]['TE'][lcntr]
                TE_der2 = cl_deriv_dic[p2]['TE'][lcntr]


            if with_lensing:
                PP_der1, TPhi_der1, EPhi_der1 = cl_deriv_dic[p]['PP'][lcntr], cl_deriv_dic[p]['Tphi'][lcntr], cl_deriv_dic[p]['Ephi'][lcntr]
                PP_der2, TPhi_der2, EPhi_der2 = cl_deriv_dic[p2]['PP'][lcntr], cl_deriv_dic[p2]['Tphi'][lcntr], cl_deriv_dic[p2]['Ephi'][lcntr]
            else:
                PP_der1 = PP_der2 = 0.
                TPhi_der1 = TPhi_der2 = 0. 
                EPhi_der1 = EPhi_der2 = 0.

            if null_TT: TT_der1 = TT_der2 = TPhi_der1 = TPhi_der2 = 0
            if null_EE: EE_der1 = EE_der2 = EPhi_der1 = EPhi_der2 = 0
            if null_TE: TE_der1 = TE_der2 = 0
            if null_PP: PP_der1 = PP_der2 = 0

            fprime1_l_vec = get_cov(TT_der1, EE_der1, TE_der1, PP_der1, TPhi_der1, EPhi_der1)
            fprime2_l_vec = get_cov(TT_der2, EE_der2, TE_der2, PP_der2, TPhi_der2, EPhi_der2)

            curr_val = np.trace( np.dot( np.dot(inv_COV_mat_l, fprime1_l_vec), np.dot(inv_COV_mat_l, fprime2_l_vec) ) )

            F[pcnt2,pcnt] += curr_val

    return F   
########################################################################################################################

def fn_fix_params(F_mat, param_names, fix_params):

    #remove parameters that must be fixed    
    F_mat_refined = []
    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 in fix_params or p2 in fix_params: continue
            F_mat_refined.append( (F_mat[pcntr2, pcntr1]) )

    totparamsafterfixing = int( np.sqrt( len(F_mat_refined) ) )
    F_mat_refined = np.asarray( F_mat_refined ).reshape( (totparamsafterfixing, totparamsafterfixing) )

    param_names_refined = []
    for p in param_names:
        if p in fix_params: continue
        param_names_refined.append(p)


    return F_mat_refined, param_names_refined

########################################################################################################################

def fn_add_prior(F_mat, param_names, prior_dic):

    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 == p2 and p1 in prior_dic:
                prior_val = prior_dic[p1]
                F_mat[pcntr2, pcntr1] += 1./prior_val**2.

    return F_mat

########################################################################################################################

def fn_get_ellipse_specs(COV, howmanysigma = 1):
    """
    Refer https://arxiv.org/pdf/0906.4123.pdf
    """
    assert COV.shape == (2,2)
    confsigma_dic = {1:2.3, 2:6.17, 3: 11.8}

    sig_x2, sig_y2 = COV[0,0], COV[1,1]
    sig_xy = COV[0,1]
    
    t1 = (sig_x2 + sig_y2)/2.
    t2 = np.sqrt( (sig_x2 - sig_y2)**2. /4. + sig_xy**2. )
    
    a2 = t1 + t2
    b2 = t1 - t2

    a = np.sqrt(a2)
    b = np.sqrt(b2)

    t1 = 2 * sig_xy
    t2 = sig_x2 - sig_y2
    theta = np.arctan2(t1,t2) / 2.
    
    alpha = np.sqrt(confsigma_dic[howmanysigma])
    
    #return (a*alpha, b*alpha, theta)
    return (a*alpha, b*alpha, theta, alpha*(sig_x2**0.5), alpha*(sig_y2**0.5))

########################################################################################################################

def fn_get_Gaussian(mean, sigma, minx, maxx, delx):

        x = np.arange(minx, maxx, delx)

        #return x, 1./(2*np.pi*sigma)**0.5 * np.exp( -(x - mean)**2. / (2 * sigma**2.)  )
        return x, np.exp( -(x - mean)**2. / (2 * sigma**2.)  )

########################################################################################################################

def fn_get_nl(els, rms_map_T, rms_map_P = None, fwhm = None, Bl = None, elknee_t = -1, alphaknee_t = 0, elknee_p = -1, alphaknee_p = 0):
    """
    compute nl - white noise + beam
    """

    if rms_map_P == None:
        rms_map_P = rms_map_T * 1.414

    if fwhm is not None:
        fwhm_radians = np.radians(fwhm/60.)
        #Bl = np.exp((-fwhm_radians**2.) * els * (els+1) /2.35)
        sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
        sigma2 = sigma ** 2
        Bl_gau = np.exp(els * (els+1) * sigma2)            

    if Bl is None:
        Bl = Bl_gau

    rms_map_T_radians = rms_map_T * np.radians(1/60.)
    rms_map_P_radians = rms_map_P * np.radians(1/60.)

    nl_TT = (rms_map_T_radians)**2. * Bl
    nl_PP = (rms_map_P_radians)**2. * Bl

    if elknee_t != -1.:
        nl_TT = np.copy(nl_TT) * (1. + (elknee_t * 1./els)**alphaknee_t )
    if elknee_p != -1.:
        nl_PP = np.copy(nl_PP) * (1. + (elknee_p * 1./els)**alphaknee_p )

    return Bl, nl_TT, nl_PP

########################################################################################################################

def fn_pad(el, cl):
    reclen_padding_zeros = max(el) - len(cl) + 1
    cl = np.concatenate( (cl, np.zeros(reclen_padding_zeros)) )
    cl[cl == 0] = max(cl) * 1e6 #some large number        
    return cl

########################################################################################################################

def fn_fisher_forecast_Aphiphi(els, Cl_deriv_dic, delta_Cl_dic, params, pspectra_to_use, min_l = 0, max_l = 6000):

    npar = len(params)
    F = np.zeros([npar,npar])
    #els = np.arange( len( delta_Cl_dic.values()[0] ) )

    pspectra_to_use_full = np.asarray( ['PP'] )

    for lcntr, l in enumerate( els ):

        if l<min_l or l>max_l:
            continue

        PP = delta_Cl_dic['PP'][lcntr]
        COV_mat_l = PP**2.
        COV_mat_l = np.mat( COV_mat_l )
        Cinv_l = sc.linalg.pinv2(COV_mat_l) #made sure that COV_mat_l * Cinv_l ~= I
        #print l, p, p2, fprime1_l_vec, fprime2_l_vec, COV_mat_l

        pspec_combinations = []
        for X in pspectra_to_use:
            for Y in pspectra_to_use:
                xind = np.where(pspectra_to_use_full == X)[0][0]
                yind = np.where(pspectra_to_use_full == Y)[0][0]
                if [Y,X, yind, xind] in pspec_combinations: continue
                pspec_combinations.append([X, Y, xind, yind])

        param_combinations = []
        for pcnt,p in enumerate(params):
            for pcnt2,p2 in enumerate(params):
                ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
                param_combinations.append([p,p2, pcnt, pcnt2])

        for (p,p2, pcnt, pcnt2) in param_combinations:
            for (X,Y, xind, yind) in pspec_combinations:

                der1 = np.asarray( [Cl_deriv_dic[p]['PP'][lcntr]] )
                der2 = np.asarray( [Cl_deriv_dic[p2]['PP'][lcntr]] )

                fprime1_l_vec = np.zeros(len(der1))
                fprime2_l_vec = np.zeros(len(der2))

                fprime1_l_vec[xind] = der1[xind]
                fprime2_l_vec[yind] = der2[yind]

                #if l > 100:
                #    from IPython import embed; embed()

                curr_val = np.dot(fprime1_l_vec, np.dot( Cinv_l, fprime2_l_vec ))

                F[pcnt2,pcnt] += curr_val

    return F    

########################################################################################################################
########################################################################################################################
########################################################################################################################
