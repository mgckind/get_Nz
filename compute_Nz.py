#!/usr/bin/python
#Code to compute N(z) from the sparse representation
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
import sys, os
import argparse


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print '\n*****************'
        sys.stderr.write('error: %s\n' % message)
        print '*****************\n'
        self.print_help()
        sys.exit(2)


parser = MyParser(description='Compute N(z) for selected galaxies from \
SVA1 Gold 1.0 catalogs by using cuts on attributes, a pre-defined mask or a list \
of coadd objects ids *Comments/questions: Matias Carrasco Kind mcarras2@illinois.edu*')
parser.add_argument("cutsfile", help="Textfile file with cuts, check default.cuts for more info")
parser.add_argument("-m", "--mask",
                    help="Numpy file with positional index mask, usually created using cutsfile, overrides cutsfile")
parser.add_argument("-i", "--ids",
                    help="Text file (ASCII) with COADD_OBJECTS_ID to be matched (slow), it returns the mask to be "
                         "used with -m option")
parser.add_argument("--path", help="Path to the catalogs (default is " + os.getcwd() + ")", default=os.getcwd())
parser.add_argument("--root", help="Root of output filenames")
parser.add_argument("--return_pdfs", action="store_true", help="Return a numpy file with the original PDFs for the "
                                                               "selected objects, the last row  "
                                                               "is the redshift positions")
parser.add_argument("--plot", action="store_true", help="Plot the N(z)")
args = parser.parse_args()

if args.mask != None and args.ids != None:
    print '\n********************'
    print 'Error: either --mask or --ids is acceptable, not both'
    print '********************\n'
    parser.print_help()
    sys.exit(2)

if args.path[-1] != '/': args.path = args.path + '/'

F = open(args.cutsfile, 'r')
FM = F.readlines()
F.close()

Mk = {}
for Line in FM:
    if Line[0] == '#' or Line[0] == '\n': continue
    line = Line.strip().split()
    name = line[0]
    vals = line[2][1:-1]
    if ',' in vals:
        Mk[name] = [float(vals.split(',')[0]), float(vals.split(',')[1])]
    else:
        Mk[name] = [float(vals)]

if args.mask == None:
    if args.ids == None:
        print '\n***  Creating mask...\n'
        for i in xrange(len(Mk)): print Mk.keys()[i], '=', Mk.values()[i]
        print
    else:
        print '\n***  Creatin mask from ids from file (slow):  ' + args.ids + ' ...\n'
else:
    print '\n***  Reading mask from : ' + args.mask + ' ...\n'

range_z = '%.1f_%.1f' % (Mk['TPZ_ZPHOT'][0], Mk['TPZ_ZPHOT'][1])

# Name of the catalog with few columns,
# This includes, RA,DEC, MAG_AUTO_I, MAG_DETMODEL_[G,R,I,Z], MODEST_CLASS, COADD_OBJECTS_ID
# and TPZ_ZPHOT (the mean of the redshift pdf) and TPZ_ZCONF (similar to ODDS in BPZ is a 0 to 1 value)
# quantifying the shape of the pdf around the mean, usually one picks value with ZCONF > 0.5
fitsfile = args.path + 'sva1_gold_1.0_short_sg.fits'

# Return the mask given the cuts defined in Mk, output_keys are extra returns
# (optional) masked, in that case those values will be: e.g.: Bmask['RA'] or
# Bmask['TPZ_ZPHOT']

if args.mask == None:
    if args.ids == None:
        Bmask = mask.get_mask(Mk, fitsfile, output_keys=['TPZ_ZPHOT'])
        masked_data = Bmask['MASK']
        maskout = 'masked_data_' + range_z
    else:
        Bmask = mask.get_mask(Mk, fitsfile, output_keys=['TPZ_ZPHOT'], input_ids=args.ids)
        masked_data = Bmask['MASK']
        maskout = 'masked_data_ids'
else:
    masked_data = load(args.mask)
    maskout = 'masked_data'

if args.root != None:
    maskout = args.root+'_'+maskout

# Save used mask
save(maskout, masked_data)

Ngal = len(masked_data)
print '*** Number of galaxies : ', Ngal
print

# Read the sparse representation file
# and get the indexes indicated by Bmask
# And also reads the header
print '*** Reading Sparse files...\n'
P2 = pf.open(args.path + 'sva1_gold_1.0.Psparse_all.fits')
zz = P2[1].data.field('redshift')
SP = P2[2].data.field('Sparse_indices')[masked_data]
P2.close()
head = ps.read_header(args.path + 'sva1_gold_1.0.Psparse_all.fits')
z = head['z']
dz = zz[1] - zz[0]

# This functions recovers the individual pdfs for the selected data
if args.return_pdfs:
    print '*** Extracting original PDFs (slow)...\n'
    P = []
    for kk in xrange(len(masked_data)):
        rep_pdf = ps.reconstruct_pdf_int(SP[kk], head)
        P.append(rep_pdf)
    P.append(z)
    if args.mask == None:
        if args.ids == None:
            pdfout = 'PDFs_' + range_z + '_tpz'
        else:
            pdfout = 'PDFs_masked_ids_tpz'
    else:
        pdfout = 'PDFs_masked_tpz'
    if args.root != None:
        pdfout = args.root + '_' + pdfout
    save(pdfout, array(P))


print '*** Creating Dictionary\n'
A = ps.create_voigt_dict(head['z'], head['mu'], head['Nmu'], head['sig'], head['Nsig'], head['Nv'])
VALS = linspace(0, 1, head['Ncoef'])
dVals = VALS[1] - VALS[0]

print '*** Computing coefficients\n'
Coef = zeros(shape(A)[1])
RR = array(map(ps.get_N, SP))
for i in xrange(Ngal): Coef[RR[i, 1, :]] += dVals * RR[i, 0, :]

N_z = dot(A, Coef)
N_z = N_z / sum(N_z) / dz

# Saving dN/dz to a file
# the first column is the center of the dz bin and the second colum tne value
# sum(N_z*dz) = 1
if args.mask == None:
    if args.ids == None:
        fileout = 'N_z_' + range_z + '_tpz.txt'
    else:
        fileout = 'N_z_masked_ids_tpz.txt'
else:
    fileout = 'N_z_masked_tpz.txt'

if args.root != None:
    fileout = args.root + '_' + fileout
print '*** Output Files: \n'
print '    N_z txt : ' + fileout
print '    Numpy mask : ' + maskout+'.npy'
if args.return_pdfs:
    print '    PDFs file : ' + pdfout

with open(fileout, 'wb') as fout:
    fout.write(b'# zmid   n(z)\n')
    ngal = '# No galaxies : %d \n' % Ngal
    fout.write(ngal)
    savetxt(fout, zip(z, N_z), fmt='%.6f')


# Ploting dN/dz
#
if args.plot:
    plt.figure()
    plt.fill_between(z, N_z, alpha=0.5, facecolor='red')
    plt.xlabel(r'$z$', fontsize=20)
    plt.ylabel(r'$dn(z)/dz$', fontsize=20)
    plt.axvline(Mk['TPZ_ZPHOT'][0], linestyle='--', color='black')
    plt.axvline(Mk['TPZ_ZPHOT'][1], linestyle='--', color='black')
    plt.show()
