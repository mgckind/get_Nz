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

### Data
For now only works with the SVA Gold 1.0 catalog. 
To get the data go the the DES Wiki (Under photo-z wiki)

### Note
Check https://github.com/mgckind/SparsePz for more information about the sparse representation

#### Contact
Matias Carrasco Kind
mcarras2@illinois.edu

