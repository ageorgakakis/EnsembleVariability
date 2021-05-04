import astropy.io.fits as fits
import numpy
import sys

'''
python script that generates mock AGN data catalogues (fits files)
in different redshift and luminosity intervals. These are to be read by 
the FORTRAN routines that perform Bayesian Model Inference and constrain
the PSD and/or BH/Mstar-scaling parameters.
'''

def openDATA(filename="DATA/data_Fig5_Paolillo+2017.dat"):
    import re

    '''
    The data from Fig. 5 of Paolillo+2017 are read
    to determine the redshfift and luminosity intervals
    of each data-point. These intervals are then used
    to select mock AGN in the same redshift and luminosity
    range. 
    '''
    
    z1=[];z2=[];lx1=[];lx2=[];s2=[];es2=[];dtmax=[];dtmin=[]
    inp = open (filename,"r")
    for line in inp.readlines():    
        column=line.split()
        p=re.compile('#')
        all=column[0]
        hash=p.match(all)   
        if(not hash and float(column[0])<4):
                z1.append(float(column[0]))
                z2.append(float(column[1]))
                lx1.append(float(column[4]))
                lx2.append(float(column[5]))
                s2.append(float(column[6]))
                es2.append(float(column[7]))
                dtmax.append(float(column[11]))
                dtmin.append(float(column[12]))

    inp.close()
    z1=numpy.array(z1);z2=numpy.array(z2)
    lx1=numpy.array(lx1);lx2=numpy.array(lx2)
    s2=numpy.array(s2);es2=numpy.array(es2)
    dtmax=numpy.array(dtmax);dtmin=numpy.array(dtmin)
    return {"ZMIN": z1, "ZMAX": z2, "LGLXMIN": lx1, "LGLXMAX": lx2, "SIG2": s2, "ESIG2": es2, "DTMAXOBS": dtmax, "DTMINOBS": dtmin}
    #return {"ZMIN": z1, "ZMAX": z2, "LGLXMIN": lx1 - numpy.log10(1.6), "LGLXMAX": lx2 - numpy.log10(1.6), "SIG2": s2, "ESIG2": es2, "DTMAXOBS": dtmax, "DTMINOBS": dtmin}
                
def select(data, ZMIN, ZMAX, LGLXMIN, LGLXMAX):

    '''
    Apply the selection function that mimics the Paollilo+17
    variability sample. 
    '''
    
    hdu = fits.open("../DATA/CDF7Ms_selection.fits")
    lgf=hdu[1].data['LGFLUX057']
    prob = hdu[1].data['PROB']
    hdu.close()

    weight = numpy.interp(data['FT'], lgf, prob)
    lgl0570 =  numpy.log10(1.507) + data['LGLX']
    
    m1 = numpy.logical_and(data['Z']>ZMIN, data['Z']<ZMAX)
    m2 = numpy.logical_and(lgl0570>LGLXMIN, lgl0570<LGLXMAX)
    m1 = numpy.logical_and(m1,m2)
    m1 = numpy.logical_and(numpy.logical_not(weight==0),m1)

    return m1, weight

#D=openDATA(filename="../DATA/data_Fig5_Paolillo+2017.dat")
D=openDATA(filename="../DATA/data_Fake.dat")
hdu = fits.open("../DATA/SIM.fits")
data=hdu[1].data
data['LGLX'] = data['LGLX'] 
columns= hdu[1].columns
hdu.close()

zlist = numpy.unique(D["ZMIN"])

for z in zlist:
    m = D['ZMIN']==z

    new={}
    for key in D.keys():
        new[key]=D[key][m]

    DMAX =  new["DTMAXOBS"]
    DMIN =  new["DTMINOBS"]
    i  = numpy.argsort(DMAX)
    for key in new.keys():
        new[key]=new[key][i]

    for i in range(new["DTMAXOBS"].size):
        print(i, new['ZMIN'][i],  new['ZMAX'][i], new['LGLXMIN'][i], new['LGLXMAX'][i], new["DTMINOBS"][i], new["DTMAXOBS"][i])
        mask, weight  = select(data, new['ZMIN'][i],  new['ZMAX'][i], new['LGLXMIN'][i], new['LGLXMAX'][i])

        prihdr = fits.Header()
        prihdr['ZMIN']= (new['ZMIN'][i])
        prihdr['ZMAX']= (new['ZMAX'][i])
        prihdr['LGLXMIN']= (new['LGLXMIN'][i])
        prihdr['LGLXMAX']= (new['LGLXMAX'][i])
        prihdr['DTMINOBS']= (new['DTMINOBS'][i])
        prihdr['DTMAXOBS']= (new['DTMAXOBS'][i])
        prihdr['SIG2']= (new['SIG2'][i])
        prihdr['ESIG2']= (new['ESIG2'][i])

        c=[]
        for col in columns.names:
            if(col in ['LGLX', 'LGM', 'Z', 'WEIGHT']):
                c.append(fits.Column(name=columns[col].name, format=columns[col].format, array=data[col][mask]))
        
        c.append(fits.Column(name='WEIGHT', format='E', array=weight[mask]))
        hdu=fits.BinTableHDU.from_columns(c, header=prihdr, name="SAMPLE")
        
        print("SAMPLE_Z1_{}_Z2_{}_DMIN_{}_DMAX_{}_LGLX_{}.fits".format(new['ZMIN'][i], new['ZMAX'][i], new['DTMINOBS'][i], new['DTMAXOBS'][i], int(new['LGLXMIN'][i]*10)/10.0))
        hdu.writeto("../DATA/SAMPLE_Z1_{}_Z2_{}_DMIN_{}_DMAX_{}_LGLX_{}.fits".format(new['ZMIN'][i], new['ZMAX'][i], new['DTMINOBS'][i], new['DTMAXOBS'][i], int(new['LGLXMIN'][i]*10)/10.0), overwrite=True)


