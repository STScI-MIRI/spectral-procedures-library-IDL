The Spectral Procedures (SP) library provides a set of tools for the analysis 
of spectral data in IDL.  The tools were designed for data from infrared 
space telescopes, such as the Low-Resolution Spectrometer (LRS) on IRAS, the 
Short-Wavelength Spectrometer (SWS) on ISO, and the Infrared Spectrograph 
(IRS) on Spitzer, although they are generic in design and could be used with 
most any data set.  They are listed in groups by function (in file 
tools_group.txt) and alphabetically (tools_list.txt).

The SP tools are designed to handle spectral data arranged in a simple
two-dimensional array.  Each row of the image is one wavelength element, with 
the columns defined as follows:
  column 0 - wavelength (assumed units are microns)
  column 1 - flux density (Jy)
  column 2 - uncertainty in flux censity (Jy)
  column 3 - segment or order identifier
Additional columns may also be present, but are not used by these tools.

The last column (order or segment identifier) was a key driver of the SP 
library, as infrared data are commonly distributed in a series of discrete 
segments.  Some typical examples:

IRAS/LRS data - Two segments, blue and red, separated at about 14 um.
ISO/SWS data - Full-scan data consist of 12 segments from four discrete
  sets of detectors.
Spitzer/IRS data:
- Short-Low and Long-Low each contain three orders obtained in two apertures
  (e.g. SL1, SL2, and the SL-bonus order, which is a 1st order spectrum
  obtained when the target in the SL2 aperture).
- SL and LL data combined have six spectral segments, which are often reduced 
  to four by combining the bonus-order data with overlapping wavelengths from 
  the other orders.
- Short-High and Long-High each contain ten orders because they are
  cross-dispersed.
- Thus spectra with SH and LH combined contain 20 orders.

Whatever the source of the data, as long as column 0 is wavelength, 
column 1 is flux density, and column 2 uses a unique integer identifier to 
distinguish each segment or order, the SP tools will handle the data.

Storing, reading, and writing SP data

SP data are commonly stored in two standard formats:  simple FITS image
files, and IPAC tables (which are stored in ASCII format).  

Spectral FITS images require no special header information, as long as the 
column definitions described above are adhered to.  They can be read in
IDL using the procedure readfits.pro, which is supplied as part of the
standard IDL Astronomy Users Library.  The procedure mrdfits.pro will also
work.

Spectral FITS files often contain in their headers the following keywords:
  NSEG - which defines the number of spectral segments;
  NSEG01 (etc.) - the number of wavelength elements in a given segment;
  COL00DEF (etc.) - defines the columns, with the following standard names:
    wavelength, flux, flux_error, order (or segment).
These are required, but if they are present, then the optional procedure
sp_rdfits.pro can be used to read them in.  (As of 2020 April 1, this
routine still has bugs and should be used with caution.)

If the data are stored in IPAC table files, the following column headers
are required:
  wavelength - assumed units are microns
  flux - actually flux density; assumed units are Jy
  error - same units as flux
  order - with unique identifiers for each spectral order or segment
SP data stored in IPAC table files should be read with rd_sptbl.pro.
It has keywords to warn it to expect optional column headers for 
error and order.

The routine wr_spfits.pro will write SP data to a spectral FITS file.
The following format should be used:

IDL>  wr_spfits,'output.fits',sp_data,len_array

where output.fits is the output file name, sp_data is the spectral data 
array, and len_array is a vector containing the lengths of the spectral
segments in sp_data.  Setting len_array = -1 tells wr_spfits to figure out 
the array for itself.  These values will be saved in the FITS header with the 
NSEG keyword.  See the procedure itself for additional keywords.

Plotting SP data

The command spplot.pro is the go-to plotting tool.  It serves basically as a 
wrapper for the standard IDL procedure plot.pro.  All the plot keywords are
available and passed along to plot.pro.  Use the keyword "/over" to overplot.
The procedure spplotbat.pro allows the user to parse a file containing a
list of spectral FITS files.  The user can type "H", "h", or "?" for help.

Basic philosophy

The SP tools are based on a series of tools of increasing complexity, with
the more complex tools making use of the simpler tools.  See the file 
tools_group.txt for a list of the tools organized by their functionality.
The simpler tools are designed to be generic, while some of the more 
complex tools are designed with more specific data sets in mind, such as 
spectra from the IRAS/LRS or the Spitzer/IRS.

The documentation for the procedures is provided in the headers of the
procedures themselves, although its completeness varies a fair bit from
one procedure to the next.  That aspect of this library should be 
considered a work in progress and will be driven to some extent by input
from other users.
