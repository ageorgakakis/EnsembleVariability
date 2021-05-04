from astropy.table import Table
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
sys.path.append('../')

import ensvar as var
import astropy.io.fits as fits
#import sys
from itertools import count


def GETMBHMSTAR(LGM, A, B):
    return A+B*(LGM-10)


def hist2d_dmh_lglx(dmh, lglx, PARAMS):
  

    import scipy

    LGLX_CUT_LO = PARAMS['LGLX_CUT_LO']
    LGLX_CUT_HI = PARAMS['LGLX_CUT_HI']
    LGM_CUT_LO = PARAMS['LGM_CUT_LO']
    LGM_CUT_HI = PARAMS['LGM_CUT_HI']
    DELTA_LGMASS= PARAMS['DELTA_LGMASS']
    DELTA_LGLX = PARAMS['DELTA_LGLX']

    NUM_LGMASS =  int ( ( LGM_CUT_HI -  LGM_CUT_LO ) / DELTA_LGMASS ) 
    NUM_LGLX =  int ( ( LGLX_CUT_HI -  LGLX_CUT_LO ) / DELTA_LGLX ) 
    xedges = numpy.arange(LGLX_CUT_LO, LGLX_CUT_HI, DELTA_LGLX)
    yedges= numpy.arange(LGM_CUT_LO, LGM_CUT_HI, DELTA_LGMASS)
    X, Y = numpy.meshgrid(xedges[0:xedges.size-1] + 0.5 * ( xedges[1:xedges.size] - xedges[0:xedges.size-1] ), yedges[0:yedges.size-1] + 0.5 * ( yedges[1:yedges.size] - yedges[0:yedges.size-1] ))
    print( "X:",len(xedges)-1)
    print("Y:",len(yedges)-1)
    levels = [1.,2,3.]

    
    #print dmh
    #print lgln
    #print(LGLX_CUT_LO,LGLX_CUT_HI,DELTA_LGLX)
 


    xp = []
    yp50 = []
    yp05 = []
    yp95 = []
    
    egdesLGLX=numpy.arange(LGLX_CUT_LO,LGLX_CUT_HI,DELTA_LGLX)
    for i in range(egdesLGLX.size-1):
        mask1 =  lglx>=egdesLGLX[i]
        mask2 =  lglx<egdesLGLX[i+1]

        #print lglx, egdesLGLX[i], egdesLGLX[i+1],
        
        mask = numpy.logical_and(mask1, mask2)

        #print mask[mask].size
        
        cur = dmh[mask]
        #print cur.size, egdesLGLX[i], egdesLGLX[i+1], dmh
        q=[0,0,0]
        if(cur.size>10):
            q = numpy.percentile(cur, q=[16,50,84])
        xp.append(0.5 * (egdesLGLX[i] +egdesLGLX[i+1] ))
        yp05.append(q[0])
        yp50.append(q[1])
        yp95.append(q[2])

    xp=numpy.asarray(xp)
    yp05=numpy.asarray(yp05)
    yp50=numpy.asarray(yp50)
    yp95=numpy.asarray(yp95)


    H, xedges, yedges = numpy.histogram2d(lglx, dmh, bins=(xedges, yedges), normed=False, weights=None)
    H = H.T # Let each row list bins with common y range.
    Z=numpy.asarray(H)
    mask = Z<=0
    Z[mask]=1e-10

    Z = Z / Z.sum()

    n = 1000
    t = numpy.linspace(0, Z.max(), n)
    integral = ((Z >= t[:, None, None]) * Z).sum(axis=(1,2))

    from scipy import interpolate
    f = interpolate.interp1d(integral, t)
    t_contours = f(numpy.array([0.95, 0.68]))

    t_contours = numpy.log10(t_contours)
    
    Z=numpy.log10(Z)
    
    #print len(xedges),len(yedges)
    #print Z.shape
    
    return {"X":X, "Y":Y, "Z": Z, "xp": xp,  "yp05": yp05,  "yp50": yp50,  "yp95": yp95, "levels":t_contours}


def read4OBS(ZMIN, ZMAX):

    Models = 'Model1'
    DMINOBS=0.25;DMAXOBS=6205
    
    filename = "DATA/cdfs7Ms_cat.fits"
    hdu = fits.open(filename)
    obs=hdu[1].data
    hdu.close()
    m = numpy.logical_and(obs['REDSHIFT']>ZMIN,obs['REDSHIFT']<ZMAX)
    m = numpy.logical_and(m, obs['FB_COUNTS']>350)
    obs=obs[m]

    hdu = fits.open("DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()
    EM = var.EmpMo(data=data, ZMIN=ZMIN, ZMAX=ZMAX, PSDfun=Models, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD)
    EM.updateSTRUCT()
    EM.sigma2(DMINOBS, DMAXOBS)
    EM.DSTRUCT['LGLX'] = EM.DSTRUCT['LGLX'] + numpy.log10(1.507)
    random = numpy.random.random_sample(len(EM.DSTRUCT['WEIGHT']))
    mask =  random<EM.DSTRUCT['WEIGHT']

    return obs, EM, mask


def resamplez():

    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    fig, ax = plt.subplots(figsize=(8, 6))
            
    obs, EM, mask = read4OBS(ZMIN=0.4, ZMAX=4.0)

    bins=numpy.arange(0,4.2,0.05)
    hobs, edges = numpy.histogram(obs['REDSHIFT'], bins=bins, range=None, normed=None, weights=None, density=True)
    mobs, edges = numpy.histogram(EM.DSTRUCT['Z'][mask], bins=bins, range=None, normed=None, weights=None, density=True)

    m = mobs>0
    ratio = numpy.zeros(len(hobs))
    ratio[m] = hobs[m]/mobs[m]

    index = EM.DSTRUCT['Z'][mask] / 0.05
    index = index.astype(int)
    weights = ratio[index]

#    bins=numpy.arange(0,4.2,0.2)
#    hist, edges = numpy.histogram(obs['REDSHIFT'], bins=bins, density=False)
#    plt.hist(EM.DSTRUCT['Z'][mask], bins=bins, weights=weights,density=True, histtype='step', color="steelblue", facecolor=None, fill=True, alpha=0.7, edgecolor="black", hatch=None, linewidth=2, label="Mock AGN")
#    plt.hist(obs['REDSHIFT'], bins=bins, density=True, histtype='step', color=None, facecolor=None, edgecolor="red", alpha=1, hatch='x', linewidth=2, label=r"7\,Ms CDFS, Luo+17")

    bins=numpy.arange(41,46,0.5)
    hist, edges = numpy.histogram(numpy.log10(obs['ABS_CORR_LX']), bins=bins, density=False)
    plt.hist(EM.DSTRUCT['LGLX'][mask], weights=weights, bins=bins, density=True, histtype='step', color="steelblue", facecolor=None, fill=True, alpha=0.7, edgecolor="black", hatch=None, linewidth=2, label="Mock AGN")
    plt.hist(numpy.log10(obs['ABS_CORR_LX']), bins=bins, density=True, histtype='step', color=None, facecolor=None, edgecolor="red", alpha=1, hatch='x', linewidth=2, label=r"7Ms CDFS, Luo+17")
    plt.ylim((0,0.8))
    plt.xlim((41.5,46.0))


    bin_centers = edges[0:-1]+0.5*(edges[1:]-edges[0:-1])    
    y = hist/float(sum(hist))/0.5
    print(y)
    print(bin_centers)
    err = y * numpy.sqrt(1./hist + 1./float(sum(hist)) )
    plt.errorbar(bin_centers,y,
    yerr = err, ecolor="red",elinewidth=2,fmt='none')
    
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    xmajorLocator = MultipleLocator(0.5)
    xmajorFormatter = FormatStrFormatter('%3.1f')
    xminorLocator =  MultipleLocator(0.5)
    ymajorLocator = MultipleLocator(0.1)
    ymajorFormatter = FormatStrFormatter('%3.1f')
    yminorLocator =  MultipleLocator(0.5)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.yaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.xaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.yaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.tick_params(labelsize=15)

    ax.set_ylabel(r'Relative fraction', fontdict=font)
    #ax.set_xlabel('Redshift', fontdict=font)
    #ax.set_xlabel(r'$\log\,L_X(\rm 0.5-7\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=font)
    ax.grid()
    ax.legend(loc="upper left", fontsize="x-large")
    #ax.legend(loc="upper right", fontsize="xx-large")
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        
    plt.show()
      

def plotHIST():

    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    fig, ax = plt.subplots(figsize=(8, 6))
    
    obs, EM, mask = read4OBS(ZMIN=0.4, ZMAX=4.0)


    #bins=numpy.arange(0,4.2,0.2)
    #hist, edges = numpy.histogram(obs['REDSHIFT'], bins=bins, density=False)
    #plt.hist(EM.DSTRUCT['Z'][mask], bins=bins, density=True, histtype='step', color="steelblue", facecolor=None, fill=True, alpha=0.7, edgecolor="black", hatch=None, linewidth=2, label="Mock AGN")
    #plt.hist(obs['REDSHIFT'], bins=bins, density=True, histtype='step', color=None, facecolor=None, edgecolor="red", alpha=1, hatch='x', linewidth=2, label=r"7\,Ms CDFS, Luo+17")

    bins=numpy.arange(41,46,0.5)
    hist, edges = numpy.histogram(numpy.log10(obs['ABS_CORR_LX']), bins=bins, density=False)
    plt.hist(EM.DSTRUCT['LGLX'][mask], bins=bins, density=True, histtype='step', color="steelblue", facecolor=None, fill=True, alpha=0.7, edgecolor="black", hatch=None, linewidth=2, label="Mock AGN")
    plt.hist(numpy.log10(obs['ABS_CORR_LX']), bins=bins, density=True, histtype='step', color=None, facecolor=None, edgecolor="red", alpha=1, hatch='x', linewidth=2, label=r"7Ms CDFS, Luo+17")
    plt.ylim((0,0.8))
    plt.xlim((41.5,46.0))


    bin_centers = edges[0:-1]+0.5*(edges[1:]-edges[0:-1])    
    y = hist/float(sum(hist))/0.5
    print(y)
    print(bin_centers)
    err = y * numpy.sqrt(1./hist + 1./float(sum(hist)) )
    plt.errorbar(bin_centers,y,
    yerr = err, ecolor="red",elinewidth=2,fmt='none')
    
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    xmajorLocator = MultipleLocator(0.5)
    xmajorFormatter = FormatStrFormatter('%3.1f')
    xminorLocator =  MultipleLocator(0.5)
    ymajorLocator = MultipleLocator(0.1)
    ymajorFormatter = FormatStrFormatter('%3.1f')
    yminorLocator =  MultipleLocator(0.5)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.yaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.xaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.yaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.tick_params(labelsize=15)

    ax.set_ylabel(r'Relative fraction', fontdict=font)
    #ax.set_xlabel('Redshift', fontdict=font)
    ax.set_xlabel(r'$\log\,L_X(\rm 0.5-7\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=font)
    ax.grid()
    ax.legend(loc="upper left", fontsize="x-large")
    #ax.legend(loc="upper right", fontsize="xx-large")
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        
    plt.show()
    
    
def plotOBS():

    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    fig, ax = plt.subplots(figsize=(8, 6))


    obs, EM, mask = read4OBS(ZMIN=0.38, ZMAX=4.0)


    DICPARAM={"LGLX_CUT_LO": 0, "LGLX_CUT_HI": 4, "LGM_CUT_LO": 42, "LGM_CUT_HI": 45.1, 'DELTA_LGMASS': 0.1, 'DELTA_LGLX': 0.1}

    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    DIC2D = hist2d_dmh_lglx(EM.DSTRUCT['LGLX'][mask], EM.DSTRUCT['Z'][mask], DICPARAM)
    ax.axis([DICPARAM['LGLX_CUT_LO'], DICPARAM['LGLX_CUT_HI'], DICPARAM['LGM_CUT_LO'], DICPARAM['LGM_CUT_HI']])
    levels = DIC2D['levels'] #[1.,2,3.]
    im = ax.imshow(DIC2D['Z'], interpolation='bilinear', origin='lower',cmap=cm.Blues, extent=(DICPARAM['LGLX_CUT_LO'], DICPARAM['LGLX_CUT_HI'], DICPARAM['LGM_CUT_LO'], DICPARAM['LGM_CUT_HI']), vmin=DIC2D['levels'][0]-0.2, vmax=DIC2D['levels'][-1]+0.5, aspect='auto')
    #contour = ax.contour(DIC2D['X'], DIC2D['Y'], DIC2D['Z'], levels=levels, origin="lower", linewidths=1.2, linestyles =('-'), colors=('black'))
    #ax.hexbin(EM.DSTRUCT['Z'][mask], EM.DSTRUCT['LGLX'][mask], gridsize=30, cmap='Blues')

    ax.scatter(obs['REDSHIFT'], numpy.log10(obs['ABS_CORR_LX']), c='red', label="7Ms CDFS, Luo+17")

    xmajorLocator = MultipleLocator(0.5)
    xmajorFormatter = FormatStrFormatter('%3.1f')
    xminorLocator =  MultipleLocator(0.5)
    ymajorLocator = MultipleLocator(0.5)
    ymajorFormatter = FormatStrFormatter('%3.1f')
    yminorLocator =  MultipleLocator(0.5)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.yaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.xaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.yaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.tick_params(labelsize=15)

    ax.set_ylabel(r'$\log\,L_X(\rm 0.5-7\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=font)
    ax.set_xlabel('Redshift', fontdict=font)
    ax.grid()
    ax.legend(loc="lower right", fontsize="large")
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()
    

  
def plotSELFUN():

    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter


    hdu = fits.open("DATA/CDF7Ms_selection.fits")
    lgf=hdu[1].data['LGFLUX057']
    prob = hdu[1].data['PROB']
    hdu.close()

    fontPanel = {'family': 'serif',
                 'color': 'black',
                 'weight': 'normal',
                'size': 24,} 

    #fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.plot(10**lgf, prob, color="k", linewidth=4.5, linestyle="--")

    xmajorLocator = MultipleLocator(0.5)
    xmajorFormatter = FormatStrFormatter('%3.1f')
    xminorLocator =  MultipleLocator(0.5)
    ymajorLocator = MultipleLocator(1)
    ymajorFormatter = FormatStrFormatter('%3.1f')
    yminorLocator =  MultipleLocator(0.5)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.yaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.xaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.yaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.tick_params(labelsize=24)

    ax.set_xlabel(r'$f_X(\rm 0.5-7\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=fontPanel)
    ax.set_ylabel(r'Probability', fontdict=fontPanel)
    #ax.set_yscale('log');
    ax.set_xscale('log')
    ax.axis([1e-16, 1e-12,0,1])
    ax.yaxis.set_ticks(numpy.arange(0,1.2,0.2))
    ax.grid()
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        
    
    plt.show()

    
def plotLGLX_BH():

    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    fig, ax = plt.subplots(figsize=(7.5, 6.5))


    Models = 'Model1'
    ZMIN=0.4
    ZMAX=4.0
    DMINOBS=0.25;DMAXOBS=6205

    filename = "DATA/SIM_AIRD15_VAR_PHOT.fits"
    hdu = fits.open(filename)
    data=hdu[1].data
    hdu.close()
    EM = var.EmpMo(data=data, ZMIN=ZMIN, ZMAX=ZMAX, PSDfun=Models, MMfun=var.GETMBHMSTAR_SAV_SIGMA, LEDDfun=var.GETLGLEDD)
    EM.updateSTRUCT()
    EM.sigma2(DMINOBS, DMAXOBS)
    EM.DSTRUCT['LGLX'] = EM.DSTRUCT['LGLX'] + numpy.log10(1.593)
    random = numpy.random.random_sample(len(EM.DSTRUCT['WEIGHT']))
    mask =  random<EM.DSTRUCT['WEIGHT']
    print(len(mask), len(mask[mask]))

    dlglx = 0.5;meant=[];lglxt=[]; sigmat=[]
    bins=numpy.arange(41,47, dlglx)
    for j, lglx in enumerate(bins):
        m1 = EM.DSTRUCT['LGLX']>=lglx
        m2 = EM.DSTRUCT['LGLX']<lglx + dlglx
        m1 = numpy.logical_and(m1,m2)

        maskall = numpy.logical_and(m1,mask)
        sumw = len(maskall[maskall])

        a = numpy.array([0,0,0])
        
        mean1=0
        sigma1 = 0
        if(sumw>0):
            #a = numpy.percentile(EM.DSTRUCT['LGMBH'][maskall], q=[16,50,84])
            #mean1 = numpy.sum(EM.DSTRUCT['SIGMA2O'][maskall])/sumw
            #mean1 = numpy.sum(EM.DSTRUCT['LGMBH'][maskall])/sumw
            #sigma1 = numpy.std(EM.DSTRUCT['LGMBH'][maskall])
            mean1 = numpy.sum(EM.DSTRUCT['LGLEDD'][maskall])/sumw
            sigma1 = numpy.std(EM.DSTRUCT['LGLEDD'][maskall])
            
        sumw = numpy.sum(EM.DSTRUCT['WEIGHT'][m1])        
        mean=0
        if(sumw>0):
            #mean = numpy.sum(EM.DSTRUCT['SIGMA2O'][m1] * EM.DSTRUCT['WEIGHT'][m1])/sumw            
            mean = numpy.sum(EM.DSTRUCT['LGMBH'][m1] * EM.DSTRUCT['WEIGHT'][m1])/sumw
        
        meant.append(mean1)
        sigmat.append(sigma1)
        #meant.append(a)
        lglxt.append(lglx+dlglx/2)
        print(lglx,lglx+dlglx,mean,mean1,len(m1[m1]), sumw)
        print(lglx+dlglx/2,a)

    #print(EM.DSTRUCT['PSDNORM'])
    #t = Table()
    #t['LGLX']=EM.DSTRUCT['LGLX'][mask]
    #t['LGM']=EM.DSTRUCT['LGM'][mask]
    #t['LGMBH']=EM.DSTRUCT['LGMBH'][mask]        
    #t['SIGMA2O']=EM.DSTRUCT['SIGMA2O'][mask]
    #t['Z']=EM.DSTRUCT['Z'][mask]
    #t['NUB']=EM.DSTRUCT['NUB'][mask]
    #t['PSDNORM']=EM.DSTRUCT['PSDNORM'][mask]    
    #t.write("filename.fits", format='fits', overwrite=True)


    #DICPARAM={"LGLX_CUT_LO": 42.0, "LGLX_CUT_HI": 45.5, "LGM_CUT_LO": 4.0, "LGM_CUT_HI": 10.5, 'DELTA_LGMASS': 0.4, 'DELTA_LGLX': 0.3}
    #DICPARAM={"LGLX_CUT_LO": 42.0, "LGLX_CUT_HI": 46.0, "LGM_CUT_LO": 5.5, "LGM_CUT_HI": 10., 'DELTA_LGMASS': 0.4, 'DELTA_LGLX': 0.3}
    DICPARAM={"LGLX_CUT_LO": 42.0, "LGLX_CUT_HI": 46.0, "LGM_CUT_LO": -4, "LGM_CUT_HI": 1.5, 'DELTA_LGMASS': 0.2, 'DELTA_LGLX': 0.3}

    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    #a = - numpy.log10(EM.DSTRUCT['NUB'] / (1+numpy.log10(EM.DSTRUCT['Z']) )) - numpy.log10(86400)
    #m = a<-4
    
    #print(EM.DSTRUCT['LGLX'][m])
    #print(EM.DSTRUCT['LGMBH'][m])
    #print(EM.DSTRUCT['Z'][m])
    #print(EM.DSTRUCT['NUB'][m])

    #print(numpy.percentile(a, q=[16,50,84]))
    #DIC2D = hist2d_dmh_lglx(a[mask], EM.DSTRUCT['LGLX'][mask], DICPARAM)                                          
    DIC2D = hist2d_dmh_lglx(EM.DSTRUCT['LGLEDD'][mask], EM.DSTRUCT['LGLX'][mask], DICPARAM)
    #DIC2D = hist2d_dmh_lglx(EM.DSTRUCT['LGMBH'][mask], EM.DSTRUCT['LGLX'][mask], DICPARAM)
    ax.axis([DICPARAM['LGLX_CUT_LO'], DICPARAM['LGLX_CUT_HI'], DICPARAM['LGM_CUT_LO'], DICPARAM['LGM_CUT_HI']])
    levels = DIC2D['levels'] #[1.,2,3.]
    im = ax.imshow(DIC2D['Z'], interpolation='bilinear', origin='lower',cmap=cm.Blues, extent=(DICPARAM['LGLX_CUT_LO'], DICPARAM['LGLX_CUT_HI'], DICPARAM['LGM_CUT_LO'], DICPARAM['LGM_CUT_HI']), vmin=DIC2D['levels'][0]-0.2, vmax=DIC2D['levels'][-1]+0.5, aspect='auto')
    contour = ax.contour(DIC2D['X'], DIC2D['Y'], DIC2D['Z'], levels=levels, origin="lower", linewidths=1.2, linestyles =('-'), colors=('black'))
    

    ax.errorbar(lglxt, meant,yerr=sigmat,  ecolor="black",elinewidth=2,fmt='ko')
    #ax.scatter(EM.DSTRUCT['LGLX'][mask], EM.DSTRUCT['LGLEDD'][mask])
    #print(EM.DSTRUCT['SIGMA2O'][mask])
    #ax.plot(DIC2D['xp'], DIC2D['yp50'], color="red", linewidth=2.0, linestyle="-")
    #ax.plot(DIC2D['xp'], DIC2D['yp05'], color="red", linewidth=2.0, linestyle="--")
    #ax.plot(DIC2D['xp'], DIC2D['yp95'], color="red", linewidth=2.0, linestyle="--")

    #ax.set_yscale('log')
    #ax.axis([41, 46,-5,0])
    #x=[42,45.5]
    #y=[0.25, 0.25]
    #ax.plot(x, numpy.log10(y), color="green", linewidth=2.0, linestyle=":", label="0.25d")
    #y=[625,625]
    #ax.plot(x, numpy.log10(y), color="purple", linewidth=2.0, linestyle=":", label="0.25d")


    
    xmajorLocator = MultipleLocator(0.5)
    xmajorFormatter = FormatStrFormatter('%3.1f')
    xminorLocator =  MultipleLocator(0.5)
    ymajorLocator = MultipleLocator(1)
    ymajorFormatter = FormatStrFormatter('%3.1f')
    yminorLocator =  MultipleLocator(0.5)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.yaxis.set_tick_params(which='major', width=1.5, size=7)
    ax.xaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.yaxis.set_tick_params(which='minor', width=1.5, size=3)
    ax.tick_params(labelsize=15)

    ax.set_xlabel(r'$\log L_X(\rm 2-10\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=font)
    ax.set_ylabel(r'$\log \lambda_{EDD}$', fontdict=font)
    #ax.set_ylabel(r'$\log M_{BH}\;\;\;(M_{\odot})$', fontdict=font)
    ax.grid()
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()

    
def readpost(filename="chains/age-post_equal_weights.dat"):

    import re
    
    p1=[] ;p2=[]; p3=[]
    inp = open (filename,"r")
    for line in inp.readlines():    
        column=line.split()
        p=re.compile('#')
        all=column[0]
        hash=p.match(all)   
        if(not hash and float(column[0])<1.5):
                p1.append(float(column[0]))
                p2.append(float(column[1]))
                p3.append(float(column[2]))
                
    inp.close()
    p1=numpy.array(p1);p2=numpy.array(p2); p3=numpy.array(p3)

    return p1, p2, p3


def openDATA(filename="../DATA/data_Fig8_Paolillo+2017.dat"):
    import re
    
    z1=[];z2=[];lx=[];lx1=[];lx2=[];s2=[];es2=[];dtmax=[];dtmin=[]
    inp = open (filename,"r")
    for line in inp.readlines():    
        column=line.split()
        p=re.compile('#')
        all=column[0]
        hash=p.match(all)
        #print(all)
        if(not hash and float(column[0])<5.5):
                z1.append(float(column[0]))
                z2.append(float(column[1]))
                lx.append(float(column[3]))
                lx1.append(float(column[4]))
                lx2.append(float(column[5]))
                s2.append(float(column[6]))
                es2.append(float(column[7]))
                dtmax.append(float(column[11]))
                dtmin.append(float(column[12]))

    inp.close()
    z1=numpy.array(z1);z2=numpy.array(z2)

    lx1=numpy.array(lx1);lx2=numpy.array(lx2); lx=numpy.asarray(lx)
    s2=numpy.array(s2);es2=numpy.array(es2)
    dtmax=numpy.array(dtmax);dtmin=numpy.array(dtmin)
    #return {"ZMIN": z1, "ZMAX": z2, "LGLXMIN": lx1 - numpy.log10(1.6), "LGLXMAX": lx2 - numpy.log10(1.6), "SIG2": s2, "ESIG2": es2, "DTMAXOBS": dtmax, "DTMINOBS": dtmin, "LGLX": lx - numpy.log10(1.6)}
    return {"ZMIN": z1, "ZMAX": z2, "LGLXMIN": lx1 - numpy.log10(1.0), "LGLXMAX": lx2 - numpy.log10(1.0), "SIG2": s2, "ESIG2": es2, "DTMAXOBS": dtmax, "DTMINOBS": dtmin, "LGLX": lx - numpy.log10(1.0)}


               
def plot4paperSin():

    hdu = fits.open("DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()


    EM=[]
    Models = ['Model1', 'Model2', 'Model3', 'Model4']
    ZMIN=[0.4]
    ZMAX=[4.0]
    for Z1, Z2 in zip(ZMIN, ZMAX):
        for model in Models:
            #EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD))
            EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SAV_SIGMA, LEDDfun=var.GETLGLEDD))
            EM[-1].updateSTRUCT()

    dlglx = 0.5
    bins=numpy.arange(41,47, dlglx)
    DTMINOBS=numpy.array([0.25]);DTMAXOBS=numpy.array([6205])

    array=numpy.ndarray([bins.size, DTMINOBS.size, len(EM)])
    array.fill(1e-10)

    for imodel in range(len(EM)):
        #print(EM[imodel].DSTRUCT['LGLX'])
        EM[imodel].DSTRUCT['LGLX'] = EM[imodel].DSTRUCT['LGLX'] + numpy.log10(1.593)
        #print(EM[imodel].DSTRUCT['LGLX'])
        #sys.exit()
        for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):           
            EM[imodel].sigma2(DMINOBS, DMAXOBS)
            indeces = numpy.digitize(EM[imodel].DSTRUCT['LGLX'], bins)
            possible = numpy.unique(indeces)
            for j in possible:
                if(j>0 and j<bins.size):
                    ci = indeces == j
                    #sumw = numpy.sum(EM[imodel].DSTRUCT['WEIGHT'][ci])
                    #array[j-1,idt,imodel] = numpy.sum(EM[imodel].DSTRUCT['SIGMA2O'][ci] * EM[imodel].DSTRUCT['WEIGHT'][ci])/sumw
                    array[j-1,idt,imodel] = numpy.average(a=EM[imodel].DSTRUCT['SIGMA2O'][ci], weights=EM[imodel].DSTRUCT['WEIGHT'][ci])
                    

    # generic plotting qunatities
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                      '#f781bf', '#a65628', '#984ea3',
                      '#999999', '#e41a1c', '#dede00']

    fontPanel = {'family': 'serif',
                 'color': 'black',
                 'weight': 'normal',
                 'size': 16,} 

    # define plot
    fig, ax = plt.subplots(len(ZMIN), len(DTMINOBS), figsize=(8,7))

    # organise the data to plot
    for iz in range(len(ZMIN)):
        for idt in range(len(DTMINOBS)):
            # select the data

            a = ax
            m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[idt], D['ZMIN']==ZMIN[iz])
            xobs = D['LGLX'][m]
            xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
            yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
            yobs=D['SIG2'][m]
            print(iz, idt, xobs, xerr)
            a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o', c="k")        

            for imodel,color in enumerate(CB_color_cycle[0:len(Models)]):
                i = imodel + iz*len(Models)
                print(array.shape, idt, i, iz)
                a.plot(bins+dlglx/2, array[:,idt,i], "-", color=color, linewidth=2.0, label=Models[imodel])
        

            a.set_yscale('log');a.set_yscale('log')
            a.axis([41.5,45.5,0.005,5])
            a.xaxis.set_ticks([42,43,44,45])
            a.grid()
            if(iz==0 and idt==0):
                a.legend()
            for axis in ['top','bottom','left','right']:
                a.spines[axis].set_linewidth(2)
                a.tick_params(labelsize=16)
                a.xaxis.set_tick_params(which='major', width=1.5, size=8)
                a.yaxis.set_tick_params(which='major', width=1.5, size=8)
                a.xaxis.set_tick_params(which='minor', width=1.5, size=5)
                a.yaxis.set_tick_params(which='minor', width=1.5, size=5)

            #a.set_title(r"$\Delta T={} - {}\,d$".format(DTMINOBS[idt], DTMAXOBS[idt]), fontdict=fontPanel)
            #pos1 = a.get_position() 
            #text = "$z={}-{}$".format(ZMIN[iz], ZMAX[iz])
            #plt.figtext(pos1.x1+0.01, 0.5*(pos1.y0+pos1.y1)+0.06, r"$z={}-{}$".format(ZMIN[iz], ZMAX[iz]), rotation=270, fontdict=fontPanel)    
            #a.xaxis.set_label_coords(1.15, -0.25)
            a.set_xlabel(r'$\log L_X(\rm 0.5-7\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=fontPanel)
            #a.yaxis.set_label_coords(-0.35, -0.1)
            a.set_ylabel(r'Normalised Excess Variance', fontdict=fontPanel)

    plt.figtext(0.18, 0.18, r"$\Delta T={} - {}\,d$".format(DTMINOBS[idt], DTMAXOBS[idt]), fontdict=fontPanel)
    plt.figtext(0.18, 0.23, r'$M_{star}-M_{BH}$: Shankar+16 ', fontdict=fontPanel)    
    #plt.figtext(0.18, 0.23, r'$M_{star}-M_{BH}$: Savorgnan+16 ', fontdict=fontPanel)    
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()




                
def plot4paper():

    hdu = fits.open("DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()


    EM=[]
    #Models = ['Model1', 'Model2', 'Model3', 'Model4']
    Models = ['Model1']
    #ZMIN=[0.4, 1.03, 1.8, 2.75]
    #ZMAX=[1.03, 1.8, 2.75, 4.0]
    ZMIN=[0.4]
    ZMAX=[1.03]
    for Z1, Z2 in zip(ZMIN, ZMAX):
        for model in Models:
            #EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD))
            EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SAV_SIGMA, LEDDfun=var.GETLGLEDD))
            EM[-1].updateSTRUCT()

    dlglx = 0.5
    bins=numpy.arange(41,47, dlglx)
    #DTMINOBS=numpy.array([0.25,0.95, 0.25,0.25]);DTMAXOBS=numpy.array([45,127, 654,6205])
    DTMINOBS=numpy.array([0.25]);DTMAXOBS=numpy.array([6205])

    array=numpy.ndarray([bins.size, DTMINOBS.size, len(EM)])
    array.fill(1e-10)


    for imodel in range(len(EM)):
        EM[imodel].DSTRUCT['LGLX'] = EM[imodel].DSTRUCT['LGLX'] + numpy.log10(1.593)
        for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):           
            EM[imodel].sigma2(DMINOBS, DMAXOBS)
            for j, lglx in enumerate(bins):
                m1 = EM[imodel].DSTRUCT['LGLX']>=lglx
                m2 = EM[imodel].DSTRUCT['LGLX']<lglx + dlglx
                m1 = numpy.logical_and(m1,m2)
                sumw = numpy.sum(EM[imodel].DSTRUCT['WEIGHT'][m1])
                if(sumw>0):
                    array[j,idt,imodel] = numpy.sum(EM[imodel].DSTRUCT['SIGMA2O'][m1] * EM[imodel].DSTRUCT['WEIGHT'][m1])/sumw
                print(lglx,lglx+dlglx,array[j,idt,imodel],len(m1[m1]), sumw)


    for imodel in range(len(EM)):
        #print(EM[imodel].DSTRUCT['LGLX'])
        #EM[imodel].DSTRUCT['LGLX'] = EM[imodel].DSTRUCT['LGLX'] + numpy.log10(1.593)
        #print(EM[imodel].DSTRUCT['LGLX'])
        #sys.exit()
        for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):           
            EM[imodel].sigma2(DMINOBS, DMAXOBS)
            indeces = numpy.digitize(EM[imodel].DSTRUCT['LGLX'], bins)
            possible = numpy.unique(indeces)
            for j in possible:
                if(j>0 and j<bins.size):
                    ci = indeces == j
                    #sumw = numpy.sum(EM[imodel].DSTRUCT['WEIGHT'][ci])
                    #array[j-1,idt,imodel] = numpy.sum(EM[imodel].DSTRUCT['SIGMA2O'][ci] * EM[imodel].DSTRUCT['WEIGHT'][ci])/sumw
                    array[j-1,idt,imodel] = numpy.average(EM[imodel].DSTRUCT['SIGMA2O'][ci], weights=EM[imodel].DSTRUCT['WEIGHT'][ci])
                    print(min(EM[imodel].DSTRUCT['LGLX'][ci]),max(EM[imodel].DSTRUCT['LGLX'][ci]), array[j-1,idt,imodel])
    sys.exit()
    # generic plotting qunatities
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                      '#f781bf', '#a65628', '#984ea3',
                      '#999999', '#e41a1c', '#dede00']

    fontPanel = {'family': 'serif',
                 'color': 'black',
                 'weight': 'normal',
                 'size': 16,} 

    # define plot
    fig, ax = plt.subplots(max(2,len(ZMIN)), max(2,DTMINOBS.size), figsize=(13,13))

    # organise the data to plot
    for iz in range(len(ZMIN)):
        for idt in range(len(DTMINOBS)):
            # select the data
            print(iz, idt, ax.shape)
            a = ax[iz,idt]
            m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[idt], D['ZMIN']==ZMIN[iz])
            xobs = D['LGLX'][m]
            xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
            yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
            yobs=D['SIG2'][m]
            print(iz, idt, xobs, xerr)
            a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o', c="k")        

            for imodel,color in enumerate(CB_color_cycle[0:len(Models)]):
                i = imodel + iz*len(Models)
                print(array.shape, idt, i, iz)
                a.plot(bins+dlglx/2, array[:,idt,i], "-", color=color, linewidth=2.0, label=Models[imodel])


            a.set_yscale('log');a.set_yscale('log')
            a.axis([41.5,45.5,0.005,5])
            a.xaxis.set_ticks([42,43,44,45])
            a.grid()
            if(iz==0 and idt==0):
                a.legend()
            for axis in ['top','bottom','left','right']:
                a.spines[axis].set_linewidth(2)
                a.tick_params(labelsize=16)
                a.xaxis.set_tick_params(which='major', width=1.5, size=8)
                a.yaxis.set_tick_params(which='major', width=1.5, size=8)
                a.xaxis.set_tick_params(which='minor', width=1.5, size=5)
                a.yaxis.set_tick_params(which='minor', width=1.5, size=5)

            if(iz==0):
                a.set_title(r"$\Delta T={} - {}\,d$".format(DTMINOBS[idt], DTMAXOBS[idt]), fontdict=fontPanel)

            if(idt==3):
                pos1 = a.get_position() 
                text = "$z={}-{}$".format(ZMIN[iz], ZMAX[iz])
                plt.figtext(pos1.x1+0.01, 0.5*(pos1.y0+pos1.y1)+0.06, r"$z={}-{}$".format(ZMIN[iz], ZMAX[iz]), rotation=270, fontdict=fontPanel)    
            if(iz==len(ZMIN)-1 and idt==1):
                a.xaxis.set_label_coords(1.15, -0.25)
                a.set_xlabel(r'$\log L_X(\rm 0.5-8\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=fontPanel)
            if(iz==1 and idt==0):
                a.yaxis.set_label_coords(-0.35, -0.1)
                a.set_ylabel(r'Normalised Excess Variance', fontdict=fontPanel)
            if(idt>0):
                a.set_yticklabels([])           
            if(iz<3):
                a.set_xticklabels([])           
        
    #plt.figtext(0.15, 0.92, r'$M_{star}-M_{BH}$: Shankar+16 ', fontdict=fontPanel)    
    plt.figtext(0.15, 0.92, r'$M_{star}-M_{BH}$: Savorgnan+16 ', fontdict=fontPanel)    
    pdf = PdfPages('test.pdf')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()





def plotFORT(fortdir=None, outplotname="test.pdf"):

    import matplotlib.gridspec as gridspec
    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from astropy.io import ascii

    NPLOT=2
    D=openDATA("../DATA/data_Fig5_Paolillo+2017.dat")
    dataBH = ascii.read("{}/fort.88".format(fortdir))  
    dataXSV = ascii.read("{}/fort.93".format(fortdir))  
    LGM=dataBH['col1']
    LGBHl=dataBH['col3']
    LGBHu=dataBH['col4']
    LGBHm=dataBH['col2']

    LGLX=dataXSV['col1']
    XSVM=dataXSV['col2']
    XSVL=dataXSV['col4']
    XSVH=dataXSV['col5']

    
    bg_color = 'white'
    fg_color = 'black'
    
    #LGM=numpy.arange(10,11.5,0.1)
    bhshankar =  7.574 + 1.946 * (LGM - 11) -0.306 * (LGM-11)**2 -0.011 * (LGM-11)**2
    bhstandard =   8.31+1.31*(LGM-11)

    fig = plt.figure(1, figsize=(13, 8))

    ax=list(range(NPLOT))
    ax[0] = plt.subplot2grid((12,12), (0,7), rowspan=4, colspan=4)
    ax[1] = plt.subplot2grid((12,12), (6,7), rowspan=4, colspan=4)


    for i in range(NPLOT):
        ax[i].grid(True, ls=":", c="grey")
    
    font = {'family': 'serif',
        'color': fg_color,
        'weight': 'normal',
        'size': 16,
        }
    fonts = {'family': 'serif',
        'color': fg_color,
        'weight': 'normal',
        'size': 14,
        }

    
    ax[1].fill_between(LGM, LGBHl, LGBHu, facecolor='red', alpha=0.5)
    ax[1].plot(LGM, LGBHm, 'b:', label="variability (median)")
    ax[1].plot(LGM, bhshankar, 'g-', label="Shankar+17")
    ax[1].plot(LGM, bhstandard, 'r--', label="Savorgnan+16")
    ax[1].legend()
    ax[1].set_xlabel('$\log\,M_{star}$',  fontdict=font)
    ax[1].set_ylabel(r'$\log\,M_{BH}$',  fontdict=font)
    ax[1].set(xlim=(9.5,12.0))
    ax[1].set(ylim=(3.5,10.5))


    ax[0].fill_between(LGLX, XSVL, XSVH, facecolor='red', alpha=0.5)
    ax[0].plot(LGLX, XSVM, 'b:')
    xobs = D['LGLX']
    xerr = [D['LGLX']-D['LGLXMIN'], D['LGLXMAX']-D['LGLX']]
    yerr= [D['ESIG2'], D['ESIG2']]
    yobs=D['SIG2']        
    ax[0].errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o')        
    ax[0].axis([42.5,45.5,0.01,2])
    ax[0].set_yscale('log')
    ax[0].set_xlabel('$\log\,L_{X}(\mathrm{0.5-7\,keV})\;\;(\mathrm{erg/s})$',  fontdict=font)
    ax[0].set_ylabel('excess variance',  fontdict=font)

    
    fig.canvas.draw()
    for i in range(NPLOT):    
        for axis in ['top','bottom','left','right']:
            ax[i].spines[axis].set_linewidth(2)
            
        ax[i].xaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].yaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].xaxis.set_tick_params(which='minor', width=1.5, size=5)
        ax[i].yaxis.set_tick_params(which='minor', width=1.5, size=5)
        labels = ax[i].get_xticklabels()
        ax[i].set_xticklabels(labels, fontdict=fonts)
        labels = ax[i].get_yticklabels()
        ax[i].set_yticklabels(labels,fontdict=fonts)
 
        
    pdf = PdfPages(outplotname)
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()




def plotFORT1(fortdir=None, outplotname="test.pdf"):

    import matplotlib.gridspec as gridspec
    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from astropy.io import ascii

    FOUR=True
    NPLOT=6
    if(FOUR):
        NPLOT=10

    dataPAR = ascii.read("{}/fort.90".format(fortdir))  

    p1 = dataPAR['col1']
    p2 = dataPAR['col2']
    p3 = dataPAR['col3']
    if(FOUR):
        p4 = dataPAR['col4']
    
    bg_color = 'white'
    fg_color = 'black'
    
    fig = plt.figure(1, figsize=(13, 8))

    ax=list(range(NPLOT))
    ax[0] = plt.subplot2grid((12,12), (0,0), rowspan=3, colspan=2)
    ax[1] = plt.subplot2grid((12,12), (3,0), rowspan=3, colspan=2)
    ax[2] = plt.subplot2grid((12,12), (6,0), rowspan=3, colspan=2)    
    ax[3] = plt.subplot2grid((12,12), (3,2), rowspan=3, colspan=2)
    ax[4] = plt.subplot2grid((12,12), (6,2), rowspan=3, colspan=2)    
    ax[5] = plt.subplot2grid((12,12), (6,4), rowspan=3, colspan=2)
    if(FOUR):
        ax[6] = plt.subplot2grid((12,12), (9,0), rowspan=3, colspan=2)
        ax[7] = plt.subplot2grid((12,12), (9,2), rowspan=3, colspan=2)
        ax[8] = plt.subplot2grid((12,12), (9,4), rowspan=3, colspan=2)
        ax[9] = plt.subplot2grid((12,12), (9,6), rowspan=3, colspan=2)

    for i in range(NPLOT):
        ax[i].grid(True, ls=":", c="grey")
    
    font = {'family': 'serif',
        'color': fg_color,
        'weight': 'normal',
        'size': 16,
        }
    fonts = {'family': 'serif',
        'color': fg_color,
        'weight': 'normal',
        'size': 14,
        }

    minimum1 = numpy.min(p1)*0.9; maximum1 = numpy.max(p1) * 1.05;
    minimum2 = numpy.min(p2)*0.9; maximum2 = numpy.max(p2) * 1.05;    
    minimum3 = numpy.min(p3)*0.9; maximum3 = numpy.max(p3) * 1.05;
    if(FOUR):
        minimum4 = numpy.min(p4)*0.9; maximum4 = numpy.max(p4) * 1.05;
    

    ax[1].hexbin(p2, p1, gridsize=20, cmap='Blues')
    ax[1].set_ylabel(r'$\alpha$',  fontdict=font)
    ax[1].set(xlim=(minimum2,maximum2))
    ax[1].set(ylim=(minimum1,maximum1))    
    labels = [item.get_text() for item in ax[1].get_xticklabels()]
    print(labels)
    labels = ['']*len(labels)
    ax[1].set_xticklabels(labels)
    
    #ax[2].scatter(p2, p3, marker="o", s=8, c="red", edgecolors="black", linewidths=0.2)
    ax[2].hexbin(p2, p3, gridsize=15, cmap='Blues')
    ax[2].set_ylabel(r'$\gamma$',  fontdict=font)
    ax[2].set_xlabel(r'$\beta$',  fontdict=font)    
    ax[2].set(xlim=(minimum2,maximum2))
    ax[2].set(ylim=(minimum3,maximum3))
    labels = [item.get_text() for item in ax[2].get_xticklabels()]
    for item in ax[2].get_yticklabels():
        print(item)
    #print(labels,ax[2].get_xticklabels() )
    ax[2].set_xticklabels(labels, rotation=90)
    
    #ax[4].scatter(p1, p3, marker="o", s=8, c="red", edgecolors="black", linewidths=0.2)
    ax[4].hexbin(p1, p3, gridsize=20, cmap='Blues')
    ax[4].set_xlabel(r'$\alpha$',  fontdict=font)
    ax[4].set(xlim=(minimum1,maximum1))
    ax[4].set(ylim=(minimum3,maximum3))    
    labels = [item.get_text() for item in ax[4].get_yticklabels()]
    labels = ['']*len(labels)
    ax[4].set_yticklabels(labels)
         
    if(FOUR):
       #ax[6].scatter(p2, p4, marker="o", s=8, c="red", edgecolors="black", linewidths=0.2)
       ax[6].hexbin(p2, p4, gridsize=20, cmap='Blues')
       ax[6].set_ylabel(r'$\delta$',  fontdict=font)
       ax[6].set_xlabel(r'$\beta$',  fontdict=font)
       ax[6].set(xlim=(minimum2,maximum2))
       ax[6].set(ylim=(minimum4,maximum4))    
       labels = [item.get_text() for item in ax[6].get_xticklabels()]
       #print(labels)
       #labels = ['']*len(labels)
       #ax[1].set_xticklabels(labels)

       #ax[7].scatter(p1, p4, marker="o", s=8, c="red", edgecolors="black", linewidths=0.2)
       ax[7].hexbin(p1, p4, gridsize=20, cmap='Blues')
       ax[7].set_xlabel(r'$\alpha$',  fontdict=font)
       ax[7].set(xlim=(minimum1,maximum1))
       ax[7].set(ylim=(minimum4,maximum4))
       labels = [item.get_text() for item in ax[7].get_yticklabels()]
       labels = ['']*len(labels)
       ax[7].set_yticklabels(labels)
       
       #ax[8].scatter(p3, p4, marker="o", s=8, c="red", edgecolors="black", linewidths=0.2)
       ax[8].hexbin(p3, p4, gridsize=20, cmap='Blues')
       ax[8].set_xlabel(r'$\gamma$',  fontdict=font)
       ax[8].set(xlim=(minimum3,maximum3))
       ax[8].set(ylim=(minimum4,maximum4))    
       labels = [item.get_text() for item in ax[8].get_yticklabels()]
       labels = ['']*len(labels)
       ax[8].set_yticklabels(labels)


     
    # HISTOGRAMES
    n, bins, patches = ax[0].hist(p2, 30,  range=(minimum2, maximum2), histtype='step', weights=None, density=1)
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
    ax[0].set(xlim=(minimum2,maximum2))
    labels = [item.get_text() for item in ax[1].get_xticklabels()]
    labels = ['']*len(labels)
    ax[0].set_xticklabels(labels)
    labels = [item.get_text() for item in ax[0].get_yticklabels()]
    labels = ['']*len(labels)
    ax[0].set_yticklabels(labels)
         
    n, bins, patches = ax[3].hist(p1, 30,  range=(minimum1, maximum1), histtype='step', weights=None, density=1,  orientation="horizontal")
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
    ax[3].set(ylim=(minimum1,maximum1))
    labels = [item.get_text() for item in ax[1].get_yticklabels()]
    labels = ['']*len(labels)
    ax[3].set_yticklabels(labels)
    labels = [item.get_text() for item in ax[3].get_xticklabels()]
    labels = ['']*len(labels)
    ax[3].set_xticklabels(labels)
    

    n, bins, patches = ax[5].hist(p3, 30,  range=(minimum3, maximum3), histtype='step', weights=None, density=1,  orientation="horizontal")
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
    ax[5].set(ylim=(minimum3,maximum3))
    labels = [item.get_text() for item in ax[2].get_yticklabels()]
    labels = ['']*len(labels)
    ax[5].set_yticklabels(labels)
    labels = [item.get_text() for item in ax[5].get_xticklabels()]
    labels = ['']*len(labels)
    ax[5].set_xticklabels(labels)

    if(FOUR):
        n, bins, patches = ax[9].hist(p4, 30,  range=(minimum4, maximum4), histtype='step', weights=None, density=1,  orientation="horizontal")
        plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
        ax[9].set(ylim=(minimum4,maximum4))
        labels = [item.get_text() for item in ax[9].get_yticklabels()]
        labels = ['']*len(labels)
        ax[9].set_yticklabels(labels)
        labels = [item.get_text() for item in ax[9].get_xticklabels()]
        labels = ['']*len(labels)
        ax[9].set_xticklabels(labels)

    fig.canvas.draw()
    for i in range(NPLOT):    
        for axis in ['top','bottom','left','right']:
            ax[i].spines[axis].set_linewidth(2)
            
        ax[i].xaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].yaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].xaxis.set_tick_params(which='minor', width=1.5, size=5)
        ax[i].yaxis.set_tick_params(which='minor', width=1.5, size=5)
        labels = ax[i].get_xticklabels()
        ax[i].set_xticklabels(labels, fontdict=fonts)
        labels = ax[i].get_yticklabels()
        ax[i].set_yticklabels(labels,fontdict=fonts)
 
    pdf = PdfPages(outplotname)
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()







def plotFIT(p1=[], p2=[], ZMIN=0.4, ZMAX=1.03, PSDfun="Model4", MMfun=None, LEDDfun=None, outplotname="test.pdf", fortdir=None):

    import matplotlib.gridspec as gridspec
    import matplotlib.cm as cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    bg_color = 'white'
    fg_color = 'black'
    
    hdu = fits.open("DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()     
    EM1 = var.EmpMo(data=data, ZMIN=ZMIN, ZMAX=ZMAX, PSDfun=PSDfun, MMfun=MMfun, LEDDfun=LEDDfun)
    
    LGM=numpy.arange(10,11.5,0.1)
    LGBHm=[];LGBHl=[];LGBHu=[]
    for t in LGM:
        BH = GETMBHMSTAR(t, p1,p2)
        BHm, BHl, BHu = numpy.percentile(BH, [50,16,84])
        LGBHm.append(BHm)
        LGBHl.append(BHl)
        LGBHu.append(BHu)
    LGBHm = numpy.asarray(LGBHm)
    LGBHu = numpy.asarray(LGBHu)
    LGBHl = numpy.asarray(LGBHl)

    bhshankar =  7.574 + 1.946 * (LGM - 11) -0.306 * (LGM-11)**2 -0.011 * (LGM-11)**2
    bhstandard =   8.31+1.31*(LGM-11)

    fig = plt.figure(1, figsize=(15, 8))

    #ax0.axis('off')
    ax=range(5)
    ax[0] = plt.subplot2grid((6,12), (0,0), rowspan=2, colspan=4)
    ax[1] = plt.subplot2grid((6,12), (2,0), rowspan=4, colspan=4)
    ax[2] = plt.subplot2grid((6,12), (2,4), rowspan=4, colspan=2)
    ax[3] = plt.subplot2grid((6,12), (2,8), rowspan=4, colspan=4)
    ax[4] = plt.subplot2grid((6,12), (0,8), rowspan=2, colspan=4)
    ax[0].grid(True, ls=":", c="grey")
    ax[1].grid(True, ls=":", c="grey")
    ax[2].grid(True, ls=":", c="grey")
    ax[3].grid(True, ls=":", c="grey")
    ax[4].grid(True, ls=":", c="grey")
    
    font = {'family': 'serif',
        'color': fg_color,
        'weight': 'normal',
        'size': 18,
        }
    #if(not (xlab is None)):
    ax[1].set_xlabel("$\log_{10}(\mathrm{Normalisation})$", fontdict=font)
    #if(not (ylab is None)):        
    ax[1].set_ylabel("slope",  fontdict=font)
  
    minimum1 = numpy.min(p1)*1.; maximum1 = numpy.max(p1) * 1.0;
    n, bins, patches = ax[0].hist(p1, 30,  range=(minimum1, maximum1), histtype='step', weights=None, density=1)
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
 
    minimum2 = numpy.min(p2)*1.0; maximum2 = numpy.max(p2) * 1.0;
    n, bins, patches = ax[2].hist(p2, 30,  range=(minimum2, maximum2), histtype='step', weights=None, density=1,  orientation="horizontal")
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
    ax[2].tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labeltop=True,
        labelbottom=False) # labels along the bottom edge are off
        #labels = [item.get_text() for item in ax[2].get_xticklabels()]
        #empty_string_labels = ['']*len(labels)
        #ax[0].set_xticklabels(empty_string_labels)
    labels = [item.get_text() for item in ax[2].get_yticklabels()]
    empty_string_labels = ['']*len(labels)
    ax[2].set_yticklabels(empty_string_labels)



    #DICPARAM={"LGLX_CUT_LO": minimum2, "LGLX_CUT_HI": maximum2, "LGM_CUT_LO": minimum1, "LGM_CUT_HI": maximum1, 'DELTA_LGMASS': 1, 'DELTA_LGLX': 0.5}
    #x = p1; y = p2    
    #DIC2DALL = hist2d_dmh_lglx(y,x, DICPARAM)
    #levelsALL = DIC2DALL['levels']

    ax[1].scatter(p1, p2, marker="o", s=10, c="red", edgecolors="black", linewidths=0.2)

    
    #im = ax[1].imshow(DIC2DALL['Z'], interpolation='bilinear', origin='lower',cmap=cm.Blues, extent=(DICPARAM['LGLX_CUT_LO'], DICPARAM['LGLX_CUT_HI'], DICPARAM['LGM_CUT_LO'], DICPARAM['LGM_CUT_HI']), vmin=DIC2DALL['levels'][0]-0., vmax=DIC2DALL['levels'][-1]+0., aspect='auto')
    #contour = ax[1].contour(DIC2DALL['X'], DIC2DALL['Y'], DIC2DALL['Z'], levels=levelsALL, origin="lower", linewidths=2, linestyles =('-'), colors=('black'))

    

    ax[3].fill_between(LGM, LGBHl, LGBHu, facecolor='red', alpha=0.5, label="variability Q90")
    ax[3].plot(LGM, LGBHm, 'b:', label="variability (median)")
    ax[3].plot(LGM, bhshankar, 'g-', label="shankar+16")
    ax[3].plot(LGM, bhstandard, 'r--', label="Markoni&Hunt+03")
    ax[3].legend()
    ax[3].set_xlabel('$\log\,M_{star}$',  fontdict=font)
    ax[3].set_ylabel(r'$\log\,M_{BH}$',  fontdict=font)

    #minimum1 = numpy.min(p1)*1.; maximum1 = numpy.max(p1) * 1.0;
    n, bins, patches = ax[4].hist(EM1.DSTRUCT['LGM'], 30, histtype='step', weights=None, density=1,)
    plt.setp(patches, color='red',  hatch="/", linestyle='-', linewidth=2)
    x1,x2,_,_ = ax[3].axis()
    ax[4].set(xlim=(x1,x2))


    fig.canvas.draw()
    for i in range(5):    
        for axis in ['top','bottom','left','right']:
            ax[i].spines[axis].set_linewidth(2)
            
        ax[i].xaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].yaxis.set_tick_params(which='major', width=1.5, size=8)
        ax[i].xaxis.set_tick_params(which='minor', width=1.5, size=5)
        ax[i].yaxis.set_tick_params(which='minor', width=1.5, size=5)
        #labels = ax[i].get_xticklabels()
        #ax[i].set_xticklabels(labels, fontdict=font)
        #labels = ax[i].get_yticklabels()
        #ax[i].set_yticklabels(labels,fontdict=font)
 
        #ax[i].xaxis.set_major_locator(MultipleLocator(1))
        #ax[i].xaxis.set_major_formatter(FormatStrFormatter('%3.1f'))
        #ax[i].xaxis.set_minor_locator(MultipleLocator(0.5))
        #ax[i].yaxis.set_major_locator(MultipleLocator(1))
        #ax[i].yaxis.set_major_formatter(FormatStrFormatter('%3.1f'))
        #ax[i].yaxis.set_minor_locator(MultipleLocator(1.0))
        
    pdf = PdfPages(outplotname)
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()



def fakeDATA(p1=10+numpy.log10(0.002), p2=1.0, ZMIN=0.4, ZMAX=1.03, PSDfun="Model4", MMfun=None, LEDDfun=None, outname="FakeData_Fig8_Paolillo+2017.dat"):

    outfilehandle = open(outname,"w")
    hdu = fits.open("../SIM_AIRD15_PHOT_SIG0p4.fits")
    data=hdu[1].data
    hdu.close()     
    EM1 = var.EmpMo(data=data, ZMIN=ZMIN, ZMAX=ZMAX, PSDfun=PSDfun, MMfun=MMfun, LEDDfun=LEDDfun)
    dlglx = 1.0
    bins=numpy.arange(41,47, dlglx)
    DTMINOBS=numpy.array([0.25,0.95, 0.25,0.25]);DTMAXOBS=numpy.array([45.,127., 654.,6205.])
    
    array=numpy.ndarray([DTMINOBS.size])
    array.fill(1e-10)
    EM1.updateSTRUCT([p1,p2])
    for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):
           
        EM1.sigma2(DMINOBS, DMAXOBS)
        #print EM1.DSTRUCT['LGMBH']
        #print EM1.DSTRUCT['LGM']
        #EM1.plotPSD(DMINOBS, DMAXOBS,[100,2001,9002])
        indeces = numpy.digitize(EM1.DSTRUCT['LGLX'], bins)
        possible = numpy.unique(indeces)
        for j in possible:
            if(j>0 and j<bins.size):
                ci = indeces == j

                array[idt] = numpy.mean(EM1.DSTRUCT['SIGMA2O'][ci])
                numpoints = (EM1.DSTRUCT['SIGMA2O'][ci]).size

                c=[]
                hdr = fits.Header()
                hdr['SIG2']=numpy.mean(EM1.DSTRUCT['SIGMA2O'][ci])
                for col in ['LGLX', 'LGM', 'Z', 'LGMBH', 'SIGMA2O', 'SIGMA2M']:
                    c.append(fits.Column(name=col, format='D', array=EM1.DSTRUCT[col][ci]))
                hdu=fits.BinTableHDU.from_columns(c, header=hdr, name="SAMPLE")
                hdu.writeto("TEST_Z1_{}_Z2_{}_DMIN_{}_DMAX_{}_LGLX_{}.fits".format(ZMIN, ZMAX, DMINOBS, DMAXOBS, int(bins[j-1]*10)/10.0), overwrite=True)

                
                error=  0.5 * array[idt]
                [fake_sigma2] = numpy.random.normal(array[idt], error, 1)
                if(numpoints>100):
                    print>>outfilehandle, ZMIN,ZMAX, numpy.mean(EM1.DSTRUCT['Z'][ci]), numpy.mean(EM1.DSTRUCT['LGLX'][ci]), bins[j-1], bins[j], fake_sigma2, error, DMINOBS, DMAXOBS, numpoints ,DMAXOBS, DMINOBS
                    #print ZMIN,ZMAX, numpy.mean(EM1.DSTRUCT['Z'][ci]), numpy.mean(EM1.DSTRUCT['LGLX'][ci]), bins[j-1], bins[j], fake_sigma2, error, DMINOBS, DMAXOBS, numpoints ,DMAXOBS, DMINOBS, numpy.mean(EM1.DSTRUCT['SIGMA2O'][ci])
    outfilehandle.close()
#    fig, ax = plt.subplots(1, DTMINOBS.size, figsize=(14,5))
#    #print ax.shape
#    for iplot, a in zip(count(), ax):
#        # select the data
#        m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[iplot], D['ZMIN']==ZMIN)
#        xobs = D['LGLX'][m]
#        xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
#        yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
#        yobs=D['SIG2'][m]
#        a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o')        
#        a.plot(bins+dlglx/2, per[0,:,iplot], "--", color="red", linewidth=2.0,label="DT={} - {}".format(DTMINOBS[iplot], DTMAXOBS[iplot]))
#        #print iplot, per[2,:,iplot], per[1,:,iplot]
#        
#        a.fill_between(bins[:]+dlglx/2, per[2,:,iplot]+1e-9, per[1,:,iplot]+1e-9, color="salmon", facecolor="none", hatch="X", edgecolor="b", linewidth=0)
#        a.set_yscale('log');a.set_yscale('log')
#        a.axis([41.5,45.5,0.01,10])
#        a.grid()
#        a.legend()
#
#    pdf = PdfPages(outplotname)
#    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
#    pdf.close()        #
#
#    plt.show()


    
def plotPSD(p1=[], p2=[], ZMIN=0.4, ZMAX=1.03, PSDfun="Model4", MMfun=None, LEDDfun=None, D={}, outplotname="test.pdf"):

    #D=openDATA("FakeData_Fig8_Paolillo+2017.dat")
    
    hdu = fits.open("DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()     
    EM1 = var.EmpMo(data=data, ZMIN=ZMIN, ZMAX=ZMAX, PSDfun=PSDfun, MMfun=MMfun, LEDDfun=LEDDfun)
    
    dlglx = 0.5
    bins=numpy.arange(41,47, dlglx)
    #DTMINOBS=numpy.array([0.25,0.95, 0.25,0.25]);DTMAXOBS=numpy.array([45.,127., 654.,6205.])
    DTMINOBS=numpy.array([0.25]);DTMAXOBS=numpy.array([6205.])
    p1=p1[0:100];p2=p2[0:100]

    array=numpy.ndarray([bins.size, p1.size,DTMINOBS.size])
    array.fill(1e-10)
    for isample, cp1,cp2 in zip(count(),p1,p2):
        EM1.updateSTRUCT([cp1,cp2])
        for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):
           
            EM1.sigma2(DMINOBS, DMAXOBS)
            #EM1.plotPSD(DMINOBS, DMAXOBS,[100,2001,9002])
            indeces = numpy.digitize(EM1.DSTRUCT['LGLX'], bins)
            possible = numpy.unique(indeces)
            for j in possible:
                if(j>0 and j<bins.size):
                    ci = indeces == j
                    array[j-1,isample,idt] = numpy.average(a=EM1.DSTRUCT['SIGMA2O'][ci], weights=EM1.DSTRUCT['WEIGHT'][ci])
                    #numpy.mean(EM1.DSTRUCT['SIGMA2O'][ci])
                    #print j, isample, idt, array[j-1,isample,idt], bins[j-1], numpy.mean(EM1.DSTRUCT['LGLX'][ci]),  DMINOBS, DMAXOBS
    per = numpy.percentile(array, q=[50,16,84], axis=1)

    
    fig, ax = plt.subplots(1, DTMINOBS.size, figsize=(14,5))
    print( DTMINOBS.size)
    if( DTMINOBS.size>1):
        for iplot, a in zip(count(), ax):
            # select the data
        
            m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[iplot], D['ZMIN']==ZMIN)
            xobs = D['LGLX'][m]
            xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
            yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
            yobs=D['SIG2'][m]
        
            a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o')        
            a.plot(bins+dlglx/2, per[0,:,iplot], "--", color="red", linewidth=2.0,label="DT={} - {}".format(DTMINOBS[iplot], DTMAXOBS[iplot]))
            #print iplot, per[2,:,iplot], per[1,:,iplot]
        
            a.fill_between(bins[:]+dlglx/2, per[2,:,iplot]+1e-9, per[1,:,iplot]+1e-9, color="salmon", facecolor="none", hatch="X", edgecolor="b", linewidth=0)
            a.set_yscale('log');a.set_yscale('log')
            a.axis([41.5,45.5,0.01,10])
            a.grid()
            a.legend()
    else:
        iplot=0
        m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[0], D['ZMIN']==ZMIN)
        xobs = D['LGLX'][m]
        xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
        yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
        yobs=D['SIG2'][m]
        
        ax.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o')        
        ax.plot(bins+dlglx/2, per[0,:,iplot], "--", color="red", linewidth=2.0,label="DT={} - {}".format(DTMINOBS[iplot], DTMAXOBS[iplot]))
        
        ax.fill_between(bins[:]+dlglx/2, per[2,:,iplot]+1e-9, per[1,:,iplot]+1e-9, color="salmon", facecolor="none", hatch="X", edgecolor="b", linewidth=0)
        ax.set_yscale('log');ax.set_yscale('log')
        ax.axis([41.5,45.5,0.01,10])
        ax.grid()
        ax.legend()



    pdf = PdfPages(outplotname)
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    pdf.close()        

    plt.show()

def PSDModel1(LGMBH):
    #nub=580/(MBH/Msolar)s-1,
    #nubxPSD(nub) = 0.02
    # PSD(nub) = A/nub * (1+nub/nub)^-1->
    # nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub)

    A = numpy.log10( 2 * 2e-2) 
    lgnub = numpy.log10(580.0) - LGMBH
    bin =0.1
    numin = (1.0+0.8) / ( 86400.0 * 6205 )
    numax = (1.0+0.8) / ( 86400.0 * 0.25 )
    
    lgnu = numpy.arange(numpy.log10(numin),numpy.log10(numax),bin)

    lgpsd = A - lgnu - numpy.log10(1 + 10**(lgnu-lgnub))

    inte = 10**A * (numpy.log(numax/numin) - numpy.log( (10**lgnub+numax) / (10**lgnub+numin) ))

    print(numin,numax,LGMBH,numpy.sum(10**lgpsd *10**lgnu * numpy.log(10)*bin)/1.3 / 0.48**(1.1-1), inte/1.3/0.48**(1.1-1))

    return lgnu, lgpsd

 
def plotPSDsimple():

    fig, ax = plt.subplots(1,1)
    fonts = {'family': 'serif',
             'color': 'black',
             'weight': 'normal',
             'size': 10,
             }             
    fontt = {'family': 'serif',
             'color': 'red',
             'weight': 'bold',
             'size': 10,
             }             


    lgnu1, lgpsd1 = PSDModel1(7.2200)
    lgnu2, lgpsd2 = PSDModel1(9.1400)

    ax.plot(10**lgnu1, 10**lgpsd1, label=r"$\log M_{BH}=7.00$")
    ax.plot(10**lgnu2, 10**lgpsd2, label=r"$\log M_{BH}=9.00$")
    ax.legend(fontsize='x-small')
    ax.set_xlabel('log(frequency) (Hz)',  fontdict=fonts)
    ax.set_ylabel('log(Amplitude)',  fontdict=fonts)
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set(ylim=(5.0,10.5))
    ax.grid()
    #ax.set_title(label="Variability\n Power Spectrum Distribution\n $PSD=f(M_{BH}, \lambda_{EDD})$", fontdict=fontt, loc='center', pad=None)
    #ax.text(0.6, 0.75, "Panel 4", transform=ax.transAxes, fontsize=10,verticalalignment='bottom', bbox=props, color="black", weight="bold")
    #pdf = PdfPages(outplotname)
    #pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    #pdf.close()        

    plt.show()


    
#fakeDATA(p1=10.0+numpy.log10(0.002), p2=1.0, ZMIN=0.4, ZMAX=1.03, PSDfun="Model4", MMfun=var.GETMBHMSTAR_SHANKAR, LEDDfun=var.GETLGLEDD, outname="FakeData_Fig8_Paolillo+2017.dat")
plotFORT(".")
plotFORT1(".")
#plotLGLX_BH()
#plotPSDsimple()
#plotOBS()
#plotHIST()

#resamplez()
sys.exit()



low=numpy.array([-1,-1.0, -3.0])
upper=numpy.array([10,3.0, 0.0])
#D=openDATA("data_Fig8_Paolillo+2017.dat")
D=openDATA("data_Fig5_Paolillo+2017.dat")
#D=openDATA("FakeData_Fig8_Paolillo+2017.dat")
#p1, p2 = readpost(filename="chains/age-post_equal_weights.dat")
#p1, p2, p3 = readpost(filename="chains/f3d-post_equal_weights.dat")
#p1 = low[0] + p1 * (upper[0] - low[0])
#p2 = low[1] + p2 * (upper[1] - low[1])
#p3 = low[2] + p3 * (upper[2] - low[2])
#print(numpy.mean(p3), numpy.quantile(p3, [.1, 0.5,0.9]))
#sys.exit()

#plotFIT(p1=p1, p2=p2, ZMIN=0.4, ZMAX=4.0,  PSDfun="Model3", MMfun=var.GETMBHMSTAR_SAV, LEDDfun=var.GETLGLEDD, outplotname="t1.pdf")
#plotPSD(p1=p1, p2=p2, ZMIN=0.4, ZMAX=4.0, PSDfun="Model3", MMfun=var.GETMBHMSTAR, LEDDfun=var.GETLGLEDD, D=D, outplotname="t2.pdf")
#plotPSD(p1=, p2=p2, ZMIN=1.03, ZMAX=1.8, PSDfun="Model4", MMfun=var.GETMBHMSTAR, LEDDfun=var.GETLGLEDD, D=D, outplotname="t3.pdf")

plot4paperSin()
#plotLGLX_BH()
#plotSELFUN()
sys.exit()

hdu = fits.open("../SIM_AIRD15_PHOT_SIG0p4.fits")
data=hdu[1].data
#ledd = numpy.log10(25/1.26e38/0.002) + data['LGLX'] - data['LGM']
#m1 = data['Z']>0.5
#m2 = data['Z']<1.0
#m3 = ledd>numpy.log10(0.25)
#m3 = data['LGLX']>44.5
#m = numpy.logical_and(m1,m2)
#m = numpy.logical_and(m,m3)
#print(data['LGLX'][m])
#print(data['LGM'][m])
#print(data['LGNH'][m])

#ledd = numpy.log10(25/1.26e38/0.002) + data['LGLX'][m] - data['LGM'][m]
#ledd = data['LGM'][m]
#print(ledd)
#print(numpy.mean(ledd), numpy.percentile(ledd, [50, 16, 84]))
#hdu.close()
#sys.exit()

ZMIN=1.03
ZMAX=1.8
EM=[];EMSH=[]
Models = ['Model1', 'Model2', 'Model3', 'Model4']
ZMIN=[0.4, 1.03]
ZMAX=[1.03, 1.8]
for Z1, Z2 in zip(ZMIN, ZMAX):
    for model in Models:
        EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD))
        EM[-1].updateSTRUCT()
        #EMSH.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD))
        #EMSH[-1].updateSTRUCT()
#print(EM1.DSTRUCT)

#LGM=numpy.arange(10,11.5,0.2)
#LGBHm=[];LGBHl=[];LGBHu=[]
#for t in LGM:
#    MSTAR =  EM[-1].DSTRUCT['LGM']
#    m = numpy.logical_and(MSTAR>t, MSTAR<t+0.1)
#    BH = EM[-1].DSTRUCT['LGMBH'][m]
#    BHm, BHl, BHu = numpy.percentile(BH, [50,16,84])
#    sigma = numpy.std(BH)
#    #print(t+0.1/2, BHm, BHl, BHu, sigma)
#    LGBHm.append(BHm)
#    LGBHl.append(BHl)
#    LGBHu.append(BHu)
#LGBHm = numpy.asarray(LGBHm)
#LGBHu = numpy.asarray(LGBHu)
#LGBHl = numpy.asarray(LGBHl)


dlglx = 0.5
bins=numpy.arange(41,47, dlglx)
DTMINOBS=numpy.array([0.25,0.95, 0.25,0.25]);DTMAXOBS=numpy.array([45,127, 654,6205])

array=numpy.ndarray([bins.size, DTMINOBS.size, len(EM)])
array.fill(1e-10)
arraysh=numpy.ndarray([bins.size, DTMINOBS.size, len(EMSH)])
arraysh.fill(1e-10)
for imodel in range(len(EM)):
    for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):           
        EM[imodel].sigma2(DMINOBS, DMAXOBS)
        indeces = numpy.digitize(EM[imodel].DSTRUCT['LGLX'], bins)
        possible = numpy.unique(indeces)
        for j in possible:
            if(j>0 and j<bins.size):
                ci = indeces == j
                array[j-1,idt,imodel] = numpy.mean(EM[imodel].DSTRUCT['SIGMA2O'][ci])
        EMSH[imodel].sigma2(DMINOBS, DMAXOBS)
        indeces = numpy.digitize(EMSH[imodel].DSTRUCT['LGLX'], bins)
        possible = numpy.unique(indeces)
        for j in possible:
            if(j>0 and j<bins.size):
                ci = indeces == j
                arraysh[j-1,idt,imodel] = numpy.mean(EMSH[imodel].DSTRUCT['SIGMA2O'][ci])

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


fig, ax = plt.subplots(len(ZMIN), DTMINOBS.size, figsize=(13,13))
#print(ax.shape)
#plt.show()
#sys.exit()
fontPanel = {'family': 'serif',
             'color': 'black',
             'weight': 'normal',
             'size': 16,} 
n1,n2 = ax.shape

for iz in range(len(ZMIN)):
    for idt in range(len(DTMINOBS)):
        # select the data
        a = ax[iz,idt]
        m = numpy.logical_and(D['DTMAXOBS']==DTMAXOBS[idt], D['ZMIN']==ZMIN[iz])
        xobs = D['LGLX'][m]
        xerr = [D['LGLX'][m]-D['LGLXMIN'][m], D['LGLXMAX'][m]-D['LGLX'][m]]
        yerr= [numpy.minimum(D['ESIG2'][m], D['SIG2'][m]*.99), D['ESIG2'][m]]
        yobs=D['SIG2'][m]
        a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o', c="k")        

        for imodel,color in enumerate(CB_color_cycle[0:len(Models)]):
            i = imodel + iz*len(Models)
            a.plot(bins+dlglx/2, array[:,idt,i], "--", color=color, linewidth=2.0)
            a.plot(bins+dlglx/2, arraysh[:,idt,i], "-", color=color, linewidth=2.0, label=Models[imodel])
        

        a.set_yscale('log');a.set_yscale('log')
        a.axis([41.5,45.5,0.005,5])
        if(iz==0):
           a.set_title(r"$\Delta T={} - {}\,d$".format(DTMINOBS[idt], DTMAXOBS[idt]), fontdict=fontPanel)
        a.grid()
        a.legend()

        for axis in ['top','bottom','left','right']:
           a.spines[axis].set_linewidth(2)
           a.tick_params(labelsize=16)
           a.xaxis.set_tick_params(which='major', width=1.5, size=8)
           a.yaxis.set_tick_params(which='major', width=1.5, size=8)
           a.xaxis.set_tick_params(which='minor', width=1.5, size=5)
           a.yaxis.set_tick_params(which='minor', width=1.5, size=5)

        if(iz==len(ZMIN)-1 and idt==1):
           #a = gca()
           a.xaxis.set_label_coords(1.15, -0.1)
           a.set_xlabel(r'$\log L_X(\rm 2-10\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=fontPanel)
        if(iz==0 and idt==0):
           a.set_ylabel(r'Normalised Excess Variance', fontdict=fontPanel)
        if(idt>0):
           a.set_yticklabels([])           
        #if(idt==0):
        #   a.set_xticklabels([])
        

        a.xaxis.set_ticks([42,43,44,45])

plt.figtext(0.15, 0.94, r'$z=0.4-1.03$; $M_{star}-M_{BH}$: Savorgnan+16 ', fontdict=fontPanel)    
#plt.title("Savornian+17 Mstar-MBH relation")
#pdf = PdfPages(outplotname)
#pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
#pdf.close()        

plt.show()
sys.exit()


bhshankar =  7.574 + 1.946 * (LGM - 11) -0.306 * (LGM-11)**2 -0.011 * (LGM-11)**3
SIGMA = 0.32 - 0.1 * (LGM - 12)
print(SIGMA)
bhstandard = 8.31+1.31*(LGM-11)

fig, axs = plt.subplots(1, 3, figsize=(10,5))
axs[0].hist(p1, bins=30)
axs[1].hist(p2, bins=30)
axs[1].set_xlabel('slope')
axs[0].set_xlabel('normalisation')
axs[2].fill_between(LGM+0.1/2, LGBHl, LGBHu, facecolor='red', alpha=0.5, label="variability Q90")
axs[2].plot(LGM, LGBHm, 'b:', label="variability (median)")
axs[2].plot(LGM, bhstandard, 'r--', label="Markoni&Hunt+03")
axs[2].plot(LGM, bhstandard-0.5, 'r--')
axs[2].plot(LGM, bhstandard+0.5, 'r--')


m =  EM[-1].DSTRUCT['LGMBH'] > 6
#axs[2].scatter(EM1.DSTRUCT['LGM'][m], EM1.DSTRUCT['LGMBH'][m], marker='o', c='black', s=1)
axs[2].plot(LGM, bhshankar, 'g-', label="shankar+16")
axs[2].plot(LGM, bhshankar-SIGMA, 'g-')
axs[2].plot(LGM, bhshankar+SIGMA, 'g-')

axs[2].legend()
axs[2].set_xlabel('$\log\,M_{star}$')
axs[2].set_ylabel(r'$\log\,M_{BH}$')
plt.show()



