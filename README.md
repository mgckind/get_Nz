get_Nz
======

Compute N(z) for DES catalogs using sparse format for photo-z PDFs for SVA 1.0 Gold
Optimize to the benchmark LSS catalogs but can be easily change for any required cuts,
it ca be used wot a pre-defined mask or from a list of coadd_object_ids


### Requirements
You need :
* **numpy**
* **scipy** (optional)
* [**pyfits**](http://www.stsci.edu/institute/software_hardware/pyfits) or **astropy**
* matplotlib (for plotting)

### Use
Make sure the fits files (from wiki) are in the same folder
Check comments on the files and wiki page for more info

To run:

    python compute_Nz.py cuts_file

The template fot the cuts file is provided

For help:

    python compute_Nz.py -h


usage: compute_Nz.py [-h] [-m MASK] [-i IDS] [--path PATH] [--root ROOT]
                     [--return_pdfs] [--plot] [--save_dict]
                     [-o OUTPUTS [OUTPUTS ...]]
                     cutsfile

Compute N(z) for selected galaxies from SVA1 Gold 1.0 catalogs by using cuts
on attributes, a pre-defined mask or a list of coadd objects ids
*Comments/questions: Matias Carrasco Kind mcarras2@illinois.edu*

positional arguments:
  cutsfile              Textfile file with cuts, check default.cuts for more
                        info

optional arguments:
  -h, --help            show this help message and exit
  -m MASK, --mask MASK  Numpy file with positional index mask, usually created
                        using cutsfile, overrides cutsfile
  -i IDS, --ids IDS     Text file (ASCII) with COADD_OBJECTS_ID to be matched
                        (slow), it returns the mask to be used with -m option
  --path PATH           Path to the catalogs (default is
                        /home/matias/Dropbox/gold/get_Nz)
  --root ROOT           Root of output filenames
  --return_pdfs         Return a numpy file with the original PDFs for the
                        selected objects, the last row is the redshift
                        positions
  --plot                Plot the N(z)
  --save_dict           Save dictionary of basis for further calculations in
                        numpy format
  -o OUTPUTS [OUTPUTS ...], --outputs OUTPUTS [OUTPUTS ...]
                        If passed, these are the names (case insentitive) of
                        the output variables ( not PDFs) for output file, ex.:
                        --outputs RA DEC TPZ_ZPHOT etc.


### Data
For now only works with the SVA Gold 1.0 catalog. 
To get the data go the the DES Wiki (Under photo-z wiki)

### Note
Check https://github.com/mgckind/SparsePz for more information about the sparse representation

#### Contact
Matias Carrasco Kind
mcarras2@illinois.edu

