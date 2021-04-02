import astropy.io.fits as fits
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import ensvar as var
import sys
from itertools import count


def selectfun():
    '''
    Determinationo of the selection function of the 
    Paolillo+17 CDFS-7Ms variability sample. This 
    is approximated as net counts in the 0.5-7keV band > 350.
    The functions uses the CDFS-7Ms exposure map, which 
    can be obtained from:

    https://personal.psu.edu/wnb3/cdfs/cdfs-chandra.html

    For a given input flux the expected number of
    net counts is estimated 

    C = F * t / ECF * EEF 

    t  : exposure map of Luo+17, assumes G=1.4
    ECF: count conversion factor 6.8211e-09/2.1062 for G=1.4 and NH=8.8e19
    EEF: is set to 0.95 and corresponds to the approximate EEF of the aperture used
    by Paolillo+17 to extract counts. 

    the conversion factor has been estimated from xspec assuming a model
    
    ========================================================================
    Model wabs<1>*powerlaw<2> Source No.: 1   Active/Off
    Model Model Component  Parameter  Unit     Value
    par  comp
    1    1   wabs       nH         10^22    8.80000E-03  +/-  0.0          
    2    2   powerlaw   PhoIndex            1.40000      +/-  0.0          
    3    2   powerlaw   norm                1.00000      +/-  0.0          
    ________________________________________________________________________

    XSPEC12>flux 0.5 7
    Model Flux    2.1061 photons 

    XSPEC12>newpar 1
    1:wabs[1]:nH:1>0.0

    ========================================================================
    Model wabs<1>*powerlaw<2> Source No.: 1   Active/Off
    Model Model Component  Parameter  Unit     Value
    par  comp
    1    1   wabs       nH         10^22    0.0          +/-  0.0          
    2    2   powerlaw   PhoIndex            1.40000      +/-  0.0          
    3    2   powerlaw   norm                1.00000      +/-  0.0          
    ________________________________________________________________________

    XSPEC12>flux  0.5 7
    Model Flux    2.1509 photons (6.8211e-09 ergs/cm^2/s) range (0.50000 - 7.0000 keV)


    Results are written to DATA/CDF7Ms_selection.fits       

    '''

    # the file CDFS-7Ms-0p5to7-bin1-02.emap.fits.gz needs to be
    # downloaded from https://personal.psu.edu/wnb3/cdfs/cdfs-chandra.html
    hdu = fits.open("DATA/CDFS-7Ms-0p5to7-bin1-02.emap.fits.gz")
    data=hdu[0].data
    hdu.close()
    data = data[data>0]
    data = numpy.log10(data)
    
    hist, bin_edges = numpy.histogram(data, numpy.arange(min(data), max(data), 0.01))
    expmedian = bin_edges[0:-1] + (bin_edges[1:] - bin_edges[0:-1])/2


    ntot = numpy.sum(hist)
    sel=[];lgf=[]
    for lgflux in numpy.arange(-19.,-10,0.01):
        counts = 10**lgflux * 10**expmedian / ( 6.8211e-09/2.1062) * 0.95
        m = counts > 350
        n1 = numpy.sum(hist[m])
        sel.append(float(n1)/float(ntot))
        lgf.append(lgflux)

        
    sel  = numpy.array(sel); lgf=numpy.array(lgf)
    c1 = fits.Column(name='LGFLUX057', array=lgf, format='E')
    c2 = fits.Column(name='PROB', array=sel, format='E')
    hdr = fits.Header()
    hdr['ECF'] =  6.8211e-09/2.1062
    hdr['EEF'] = 0.95
    hdr['MINCNTS'] = 350
    hdr['COMMENT'] = "Paolillo+17 sel. function, net counts > 350"
    t = fits.BinTableHDU.from_columns([c1, c2], header=hdr)
    t.writeto('DATA/CDF7Ms_selection.fits')
    

def readDATA(filename="data_Fig8_Paolillo+2017.dat"):

    '''
    This function is reading the Paolillo+2017 excess variance vs Lx
    data points. 
    '''

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
    return {"ZMIN": z1, "ZMAX": z2, "LGLXMIN": lx1 - numpy.log10(1.0), "LGLXMAX": lx2 - numpy.log10(1.0), "SIG2": s2, "ESIG2": es2, "DTMAXOBS": dtmax, "DTMINOBS": dtmin, "LGLX": lx - numpy.log10(1.0)}


def plotmodel():

    '''
    this function is plots the excess variance vs Lx data points of Paolillo+2017
    (their Fig. 5) and compares with the model prediction. 
    '''

    
    hdu = fits.open("../DATA/SIM_AIRD15_VAR_PHOT.fits")
    data=hdu[1].data
    hdu.close()

    D=readDATA("../DATA/data_Fig5_Paolillo+2017.dat")
    
    EM=[]
    Models = ['Model1', 'Model2', 'Model3', 'Model4']
    ZMIN=[0.4]
    ZMAX=[4.0]
    for Z1, Z2 in zip(ZMIN, ZMAX):
        for model in Models:
            EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SHANKAR_SIGMA, LEDDfun=var.GETLGLEDD))
            #EM.append(var.EmpMo(data=data, ZMIN=Z1, ZMAX=Z2, PSDfun=model, MMfun=var.GETMBHMSTAR_SAV_SIGMA, LEDDfun=var.GETLGLEDD))
            EM[-1].updateSTRUCT()

    dlglx = 0.5
    bins=numpy.arange(41,47, dlglx)
    DTMINOBS=numpy.array([0.25]);DTMAXOBS=numpy.array([6205])

    array=numpy.ndarray([bins.size, DTMINOBS.size, len(EM)])
    array.fill(1e-10)

    for imodel in range(len(EM)):

        EM[imodel].DSTRUCT['LGLX'] = EM[imodel].DSTRUCT['LGLX'] + numpy.log10(1.593)

        for idt, DMINOBS, DMAXOBS in zip(count(),DTMINOBS,DTMAXOBS):           
            EM[imodel].sigma2(DMINOBS, DMAXOBS)
            indeces = numpy.digitize(EM[imodel].DSTRUCT['LGLX'], bins)
            possible = numpy.unique(indeces)
            for j in possible:
                if(j>0 and j<bins.size):
                    ci = indeces == j

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
            #print(iz, idt, xobs, xerr)
            a.errorbar(xobs, yobs, xerr=xerr, yerr=yerr, fmt='o', c="k")        

            for imodel,color in enumerate(CB_color_cycle[0:len(Models)]):
                i = imodel + iz*len(Models)
                #print(array.shape, idt, i, iz)
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

            a.set_xlabel(r'$\log L_X(\rm 0.5-8\,keV)\;\;\; (erg\,s^{-1}$)', fontdict=fontPanel)
            a.set_ylabel(r'Normalised Excess Variance', fontdict=fontPanel)

    plt.figtext(0.18, 0.18, r"$\Delta T={} - {}\,d$".format(DTMINOBS[idt], DTMAXOBS[idt]), fontdict=fontPanel)
    plt.figtext(0.18, 0.23, r'$M_{star}-M_{BH}$: Shankar+16 ', fontdict=fontPanel)    
    #plt.figtext(0.18, 0.23, r'$M_{star}-M_{BH}$: Savorgnan+16 ', fontdict=fontPanel)    

    #pdf = PdfPages('test.pdf')
    #pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    #pdf.close()        

    plt.show()


