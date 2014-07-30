# Use for masking objects from catalog
#
__author__ = "Matias Carrasco Kind"
from numpy import *
import pyfits as pf
import sys, os


def get_mask(cuts, filename, output_keys='', input_ids=''):
    if cuts.has_key('MODEST_CLASS') and cuts.has_key('TPZ_SG_CLASS'):
        print '***********************'
        print '*        ERROR        *'
        print '***********************'
        print '* Both keys '
        print '* MODEST_CLASS'
        print '* TPZ_SG_CLASS'
        print '* for S/G separation are present, you should pick only one'
        print '* and comment out the other one'

        sys.exit(0)

    G = pf.open(filename)
    print '*** Info about the file...\n'
    print G.info()
    #print G[1].columns
    if input_ids == '':
        TB = G[1].data
        mask = ones(G[1].header['NAXIS2'], dtype='bool')
        G.close()
        for k in cuts.keys():
            if k.find('-') == -1:
                if len(cuts[k]) == 2:
                    mtemp = TB.field(k) >= cuts[k][0]
                    mask *= mtemp
                    mtemp = TB.field(k) <= cuts[k][1]
                    mask *= mtemp
                else:
                    mtemp = TB.field(k) == cuts[k][0]
                    mask *= mtemp
            else:
                minus = k.find('-')
                mtemp = (TB.field(k[:minus]) - TB.field(k[minus + 1:])) > cuts[k][0]
                mask *= mtemp
                minus = k.find('-')
                mtemp = (TB.field(k[:minus]) - TB.field(k[minus + 1:])) < cuts[k][1]
                mask *= mtemp
        output = {}
        output['MASK'] = where(mask * 1. == 1)[0]
        if output_keys != '':
            for k in output_keys:
                if k.find('-') == -1:
                    output[k] = TB.field(k)[mask]
                else:
                    minus = k.find('-')
                    output[k] = TB.field(k[:minus])[mask] - TB.field(k[minus + 1:])[mask]
        del mask, TB
        return output
    else:
        TB = G[1].data
        mask = zeros(G[1].header['NAXIS2'], dtype='bool')
        G.close()
        in_ids = loadtxt(input_ids, dtype='int')
        for k in in_ids:
            mtemp = TB.field('COADD_OBJECTS_ID') == k
            mask += mtemp
        output = {}
        output['MASK'] = where(mask * 1. == 1)[0]
        if output_keys != '':
            for k in output_keys:
                if k.find('-') == -1:
                    output[k] = TB.field(k)[mask]
                else:
                    minus = k.find('-')
                    output[k] = TB.field(k[:minus])[mask] - TB.field(k[minus + 1:])[mask]
        del mask, TB
        return output


