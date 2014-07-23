# Code to compute N(z) from the sparse representation
# presented in http://adsabs.harvard.edu/abs/2014MNRAS.441.3550C
# using TPZ (http://lcdm.astro.illinois.edu/code/mlz.html)
# Comments/questions: Matias Carrasco Kind mcarras2@illinois.edu
# 
__author__ = 'Matias Carrasco Kind'
from numpy import *
import matplotlib.pyplot as plt
import pyfits as pf
import pdf_storage as ps
import mask
import sys,os


# This is the mask used to get N(z) and galaxies in general
# for sva1_gold_1.0 all the objects
# What's there is what's included
# For a range of values used a [min_val, max_val]
# For equalities used [value]
# You should use only one of the S/G keys AND COMMENT THE OTHER ONE
# For TPZ_SG_CLASS the range goes from 0 (galaxies) to 1 (stars)


#The default values are for the LSS benchmark cuts
Mk={}
Mk['MAG_AUTO_I']=[18,22.5]
Mk['RA']=[60,95]
Mk['DEC']=[-62,-40]
Mk['MODEST_CLASS']=[1] # USE MODEST_CLASS or TPZ_SG_CLASS between 0 (galaxies) and 1 (stars)
#Mk['TPZ_SG_CLASS']=[0.0,0.8] #TPZ PROB S/G Classification
Mk['TPZ_ZCONF']=[0.6,1.1]  #ZCONF, goodness of the PDF, quality cut
Mk['TPZ_ZPHOT']=[0.2,0.4] #REDSHIFT BIN FROM TPZ
Mk['MAG_DETMODEL_G-MAG_DETMODEL_R']=[0,3]
Mk['MAG_DETMODEL_R-MAG_DETMODEL_I']=[0,2]
Mk['MAG_DETMODEL_I-MAG_DETMODEL_Z']=[0,3]


print '***  Creating mask...\n'
for i in xrange(len(Mk)):print Mk.keys()[i],'=',Mk.values()[i]
print
    

range_z='%.1f_%.1f' % (Mk['TPZ_ZPHOT'][0],Mk['TPZ_ZPHOT'][1])

# Name of the catalog with few columns,
# This includes, RA,DEC, MAG_AUTO_I, MAG_DETMODEL_[G,R,I,Z], MODEST_CLASS, COADD_OBJECTS_ID
# and TPZ_ZPHOT (the mean of the redshift pdf) and TPZ_ZCONF (similar to ODDS in BPZ is a 0 to 1 value)
# quantifying the shape of the pdf around the mean, usually one picks value with ZCONF > 0.5
fitsfile='sva1_gold_1.0_short_sg.last.fits'

# Return the mask given the cuts defined in Mk, output_keys are extra returns
# (optional) masked, in that case those values will be: e.g.: Bmask['RA'] or
# Bmask['TPZ_ZPHOT']


Bmask=mask.get_mask(Mk, fitsfile, output_keys=['TPZ_ZPHOT'])
masked_data=Bmask['MASK']
save('masked_data_'+str(Mk['TPZ_ZPHOT'][0])+'_'+str(Mk['TPZ_ZPHOT'][1]),masked_data)
# You can load the mask instead, comment the previous lines and uncomment the following
# 
#masked_data=load('masked_data_0.2_0.4.npy')

Ngal=len(masked_data)
print '*** Number of galaxies : ', Ngal ; print

# Read the sparse representation file
# and get the indexes indicated by Bmask
# And also reads the header
print '*** Reading Sparse files...\n'
P2=pf.open('sva1_gold_1.0.Psparse_all.last.fits')
zz=P2[1].data.field('redshift')
SP = P2[2].data.field('Sparse_indices')[masked_data]
P2.close()
head = ps.read_header('sva1_gold_1.0.Psparse_all.last.fits')
z = head['z']
dz=zz[1]-zz[0]

# This functions recovers the individual pdfs itself
#kk=111
#rep_pdf = ps.reconstruct_pdf_int(SP[kk], head)
#plt.plot(z,rep_pdf)

print '*** Creating Dictionary\n'
A=ps.create_voigt_dict(head['z'],head['mu'],head['Nmu'],head['sig'],head['Nsig'],head['Nv'])
VALS = linspace(0, 1, head['Ncoef'])
dVals = VALS[1] - VALS[0]

print '*** Computing coefficients\n'
Coef=zeros(shape(A)[1])
RR=array(map(ps.get_N,SP))
for i in xrange(Ngal): Coef[RR[i,1,:]]+=dVals*RR[i,0,:]

N_z=dot(A,Coef)
N_z=N_z/sum(N_z)/dz

# Ploting dN/dz
#
plt.figure()
plt.fill_between(z,N_z,alpha=0.5,facecolor='red')
plt.xlabel(r'$z$',fontsize=20)
plt.ylabel(r'$dn(z)/dz$',fontsize=20)
plt.axvline(Mk['TPZ_ZPHOT'][0],linestyle='--',color='black')
plt.axvline(Mk['TPZ_ZPHOT'][1],linestyle='--',color='black')
plt.savefig('N_z_'+range_z+'_tpz.png',bbox_inches='tight')
# Saving dN/dz to a file
# the first column is the center of the dz bin and the second colum tne value
# sum(N_z*dz) = 1
with open('N_z_'+range_z+'_tpz.txt','wb') as fout:
    fout.write(b'# zmid   n(z)\n')
    ngal='# No galaxies : %d \n' %Ngal
    fout.write(ngal)
    savetxt(fout,zip(z,N_z),fmt='%.6f')
       
plt.show()
