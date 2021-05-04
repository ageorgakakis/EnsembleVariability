import astropy.io.fits as fits
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

def GETLGLLBOL(LGLX):
    return numpy.log10(25.0) + LGLX 

def GETLGLBOL_DURAS(LGLBOL):
    '''
    bolometric correction from 
    https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..73D/abstract
    '''
    a = 10.96;b=11.93;c=17.79
     
    K  = a * (1+ (LGLBOL/b)**c)

    LGLX = LGLBOL - numpy.log10(K) +  numpy.log10(3.826e33)
    LGLBOLCGS = LGLBOL +  numpy.log10(3.826e33)

    return LGLX, LGLBOLCGS, K

def GETLGLEDD(LGMBH, LGLX):

    '''
    estimate the bolometric correction
    and then the Eddington ratio. 
    In the case of the Duras+21 KBOL
    a look-up table is created and Lbol
    values are interpolated from that. 

    '''
    # comment out for constand KBOL
    #LGLBOL = GETLGLLBOL(LGLX)

    # the following three lines are for the Duras+21
    # KBOL
    LGLBOL = numpy.arange(5,16,0.1)
    LGLX1, LGLBOLCGS1, K = GETLGLBOL_DURAS(LGLBOL)            
    LGLBOL =  numpy.interp(LGLX, LGLX1, LGLBOLCGS1)

    return LGLBOL - numpy.log10(1.26)- 38.0 - LGMBH

    
def GETMBHMSTAR(LGM, Z, PARAM):
    A = PARAM[0]
    B = PARAM[1]
    return A+B*(LGM-10)

def GETMBHMSTAR_SIGMA(LGM, Z, PARAM):
    SIGMA = 0.0*LGM+0.5
    dev = numpy.random.normal(loc=0, scale=SIGMA)
    A = PARAM[0]
    B = PARAM[1]
    return A+B*(LGM-10)+dev


def GETMBHMSTAR_SAV(LGM, Z, *PARAM):
    return 8.35+1.31*(LGM-11)

def GETMBHMSTAR_SAV_SIGMA(LGM, Z, *PARAM):
    MBH = GETMBHMSTAR_SAV(LGM, Z, *PARAM)
    SIGMA = 0.0*LGM+0.5
    dev = numpy.random.normal(loc=0, scale=SIGMA)
    MBH = MBH + dev
    return MBH


def GETMBHMSTAR_SHANKAR(LGM, Z, *PARAM):
    A = 7.574
    B = 1.946
    C = -0.306
    D = -0.011
    Y = LGM-11.0
    Y2 = Y * Y
    Y3 = Y2 * Y
    return A + B*Y + C*Y2 + D*Y3


def GETMBHMSTAR_SHANKAR_SIGMA(LGM, Z, *PARAM):
    MBH  = GETMBHMSTAR_SHANKAR(LGM, Z)
    SIGMA = 0.32 - 0.1 * (LGM - 12)
    dev = numpy.random.normal(loc=0, scale=SIGMA)
    MBH = MBH + dev    
    return MBH


class EmpMo(object):

    def __init__(self, data=None, ZMIN=0, ZMAX=1.0, PSDfun=None, MMfun=None, LEDDfun=None):

        assert not(PSDfun is None)
        assert not(MMfun is None)
        assert not(LEDDfun is None)
        assert not(data is None)
        
        assert PSDfun in ["Model1", "Model2", "Model3", "Model4"]        
        self.PSDmethod = {'Model1': self.__Model1, 'Model2': self.__Model2, 'Model3': self.__Model3, 'Model4': self.__Model4}
         
        self.ZMIN = ZMIN
        self.ZMAX = ZMAX
        self.DSTRUCT=None
        self.__select(data)
        self.PSDfun = PSDfun
        self.MMfun = MMfun
        self.LEDDfun = LEDDfun

        self.biasC = 1.3
        self.biasBETA = 1.1

    def __select(self,data):

        
        # the lines below can be uncommented to take into 
        # account the Paolillo+17 CDFS-7Ms selection function
        # They create weights for each mock AGN which
        # represent the probability of beinh inlcuded in the
        # Paolillo+17 CDFS-7Ms variability sample. 
        #if (not os.path.isfile('DATA/CDF7Ms_selection.fits')):
        #    self.__selectfun()
        #
        #hdu = fits.open("DATA/CDF7Ms_selection.fits")
        #lgf=hdu[1].data['LGFLUX057']
        #prob = hdu[1].data['PROB']
        #hdu.close()
        #weight = numpy.interp(data['FT'], lgf, prob)

        # this is using weight = 1, i.e. no selection function
        weight = numpy.ones(len(data['Z']))
        
        #print(weight)
        #print(self.ZMIN, self.ZMAX)
        
        m = numpy.logical_and(data['Z']>self.ZMIN, data['Z']<self.ZMAX)
        m = numpy.logical_and(numpy.logical_not(weight==0),m)

        data1 = data[m]
        weight = weight[m]
        if(self.DSTRUCT is None):
            self.DSTRUCT = {"Z": data1['Z'], "WEIGHT": weight, "LGM": data1['LGM'], "LGLX": data1['LGLX'], "LGMBH": numpy.empty(data1['LGLX'].size), "LGLEDD": numpy.empty(data1['LGLX'].size), "NUB": numpy.empty(data1['LGLX'].size), "PSDNORM": numpy.empty(data1['LGLX'].size), "SIGMA2O": numpy.empty(data1['LGLX'].size), "SIGMA2M": numpy.empty(data1['LGLX'].size)}
        else:
            TMP = {"Z": data1['Z'], "LGM": data1['LGM'], "LGLX": data['LGLX']}
            self.DSTRUCT.update(TMP)

 
    def updateSTRUCT(self,MMPARAMS=None):
        self.__getMBH(MMPARAMS)
        self.__getLGLEDD()
        self.PSDmethod[self.PSDfun]()
        
    def __getMBH(self, *MMPARAMS):
        assert self.DSTRUCT['LGM'].size>0 and self.DSTRUCT['Z'].size>0
        TMP={'LGMBH': self.MMfun(self.DSTRUCT['LGM'], self.DSTRUCT['Z'],*MMPARAMS)}
        self.DSTRUCT.update(TMP)
        
    def __getLGLEDD(self):
        assert self.DSTRUCT['LGMBH'].size>0 and self.DSTRUCT['LGLX'].size>0
        TMP={'LGLEDD': self.LEDDfun(self.DSTRUCT['LGMBH'], self.DSTRUCT['LGLX'])}
        self.DSTRUCT.update(TMP)
        
    def __getLGLBOL(self):
        assert self.DSTRUCT['LGLX'].size>0
        return GETLGLLBOL(self.DSTRUCT['LGLX'])
        

    def __Model1(self):
        '''
        PSD Model 1 of Paolillo+17                
        nub=580/(MBH/Msolar)s-1,
        nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub) = 0.02
        '''
        if(self.DSTRUCT['LGMBH'].size==0 or self.DSTRUCT['LGLEDD'].size==0):
            self.__getMM()
        A = numpy.log10( 2e-2) + 0.0*self.DSTRUCT['LGMBH']
        nub = numpy.log10(580.0) - self.DSTRUCT['LGMBH']
        TMP={"NUB": 10**nub, "PSDNORM": 10**A}
        self.DSTRUCT.update(TMP)
        
    def __Model2(self):
        '''
        PSD Model 2 of Paolillo+17                
        nub=580/(MBH/Msolar)s-1,

        
        #nub=580/(MBH/Msolar)s-1,
        #nubxPSD(nub) = 0.02
        # PSD(nub) = A/nub * (1+nub/nub)^-1->
        # nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub)
        '''

        if(self.DSTRUCT['LGMBH'].size==0 or self.DSTRUCT['LGLEDD'].size==0):
            self.__getMM()
        A = numpy.log10(0.02)
        LGLBOL = self.DSTRUCT['LGLEDD'] + numpy.log10(1.26)+ 38.0 + self.DSTRUCT['LGMBH']
        nub = numpy.log10(200./86400.) + (LGLBOL-44.0) - 2.0*(self.DSTRUCT['LGMBH']-6.0)
        TMP={"NUB": 10**nub, "PSDNORM": 10**A}

        self.DSTRUCT.update(TMP)
        
        
        
    def __Model3(self):
        #nub=580/(MBH/Msolar)s-1,
        #nubxPSD(nub) = 3e-3 x LEdd^-0.8
        # PSD(nub) = A/nub * (1+nub/nub)^-1->
        # nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub)
        if(self.DSTRUCT['LGMBH'].size==0 or self.DSTRUCT['LGLEDD'].size==0):
            self.__getMM()
        A = numpy.log10(3e-3)-0.8*self.DSTRUCT['LGLEDD']
        nub = numpy.log10(580) - self.DSTRUCT['LGMBH']
        TMP={"NUB": 10**nub, "PSDNORM": 10**A}
        self.DSTRUCT.update(TMP)

    def __Model4(self):
        #nub=580/(MBH/Msolar)s-1,
        #nubxPSD(nub) = 3e-3 x LEdd^-0.8
        # PSD(nub) = A/nub * (1+nub/nub)^-1->
        # nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub)
        if(self.DSTRUCT['LGMBH'].size==0 or self.DSTRUCT['LGLEDD'].size==0):
            self.__getMM()
        A = numpy.log10(3e-3)-0.8*self.DSTRUCT['LGLEDD']
        LGLBOL = self.DSTRUCT['LGLEDD'] + numpy.log10(1.26)+ 38.0 + self.DSTRUCT['LGMBH']
        nub = numpy.log10(200./86400.) + (LGLBOL-44.0) - 2.0*(self.DSTRUCT['LGMBH']-6.0)
        TMP={"NUB": 10**nub, "PSDNORM": 10**A}
        self.DSTRUCT.update(TMP)


    def __numin(self,DMAXOBS):
        return  (1.0+self.DSTRUCT['Z']) / ( 86400.0 * DMAXOBS );

    def __numax(self,DMINOBS):
        return (1.0+self.DSTRUCT['Z']) / ( 86400.0 * DMINOBS );

    def plotPSD(self, DMINOBS, DMAXOBS, isample=[1000,1001]):

        numin1 = self.__numin(DMAXOBS)
        numax1 = self.__numax(DMINOBS)        
        lgnu=numpy.arange(-9,-1,0.1)

        for i in isample:
            psd = self.DSTRUCT['PSDNORM'][i] * 1./10**lgnu * 1./(1+10**lgnu / self.DSTRUCT['NUB'][i])                

            plt.plot(lgnu, numpy.log10(psd), linewidth=2.0, label=r"$\log MBH={}, z={}, \sigma^2={}, \log\lambda={}$".format(int(self.DSTRUCT['LGMBH'][i]*100)/100.,int(self.DSTRUCT['Z'][i]*100)/100.,int(self.DSTRUCT['SIGMA2O'][i]*100)/100 ,int(self.DSTRUCT['LGLEDD'][i]*100)/100.))
            _,_,ymin, ymax = plt.axis()
            y=numpy.array([ymin,ymax])
            xmin=numpy.log10(0*y+numin1[i])
            xmax=numpy.log10(0*y+numax1[i])
            xnu=numpy.log10(0*y+self.DSTRUCT['NUB'][i])
            plt.plot(xmin,y, "-", color="red", linewidth=1.0)
            plt.plot(xmax,y, "-", color="red", linewidth=1.0)
            plt.plot(xnu,y, "-", color="blue", linewidth=1.0)
            
            plt.title("{}, DMIN:{}, DMAX:{} days".format(self.PSDfun, DMINOBS, DMAXOBS))
        plt.xlabel(r"$\log \nu$")
        plt.ylabel(r"PSD")
        plt.legend()
        plt.show()
        
    def sigma2(self, DMINOBS, DMAXOBS):
        numin1 = self.__numin(DMAXOBS)
        numax1 = self.__numax(DMINOBS)        
        

        sigma2_model = self.DSTRUCT['PSDNORM'] * (numpy.log(numax1/numin1) - numpy.log( (self.DSTRUCT['NUB']+numax1) / (self.DSTRUCT['NUB']+numin1) ) )

        TMP={'SIGMA2O': sigma2_model / self.biasC / 0.48**(self.biasBETA-1), "SIGMA2M": sigma2_model}


        self.DSTRUCT.update(TMP)


        

        
