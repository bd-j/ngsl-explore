 STIS Next Generation Spectral Library 
 (Version 1, March 2008) 
 
Don Lindler, Sigma Space Corporation, Don.J.Lindler@nasa.gov 
Sara R Heap, NASA/Goddard Space Flight Center, Sara.R.Heap@nasa.gov 
 
The STIS Next Generation Spectral Library contains STIS spectra of 370 stars observed 
in HST programs 9088, 9786, and 10222.  Each spectrum includes spectral segments 
from gratings G230LB, G430L, and G750L merged to form a single spectrum covering 
~0.2-1.0 µ. NGSL spectra were constructed via custom pipeline processing (section 2) 
followed by corrections for in-order scatter by grating G230LB  (section 3.1) and 
corrections for throughput variations caused by the mis-centering of the target in the 0.2 
arcsecond-wide slit. (section 3.2). A summary of the observations and derived stellar 
parameters is given in section 4. 
 
1.  Fil e Content and Format 
 
Each data file is a FITS binary table.  The primary FITS header contains the following 
information: 
 
FITS File primary header 
 
 
 
 
 
SIMPLE  =                    T /image conforms to FITS standard 
BITPIX  =                   16 /bits per data value 
NAXIS   =                    0 /number of axes 
EXTEND  =                    T /file may contain extensions 
TARGNAME= 'BD+112998'          /target name 
SPECTYPE= 'composite'          /Spectra type 
TELESCOP= 'HST     '           /observatory 
INSTRUME= 'STIS    '           /instrument 
FILENAME= 'h_stis_ngsl_bd+112998_v1.fits' / 
RA      =        3713.54743480 /Right Ascension (deg) 
DEC     =        10.9976309090 /Declination (deg) 
EQUINOX =              2000.00 /equinox of celestial coord. system 
APERTURE= '52X0.2          '   / 
GRATING = 'G230LB, G430L, G750L' / 
OBSDATE = '2004-04-10'         /date of observation YYYY-MM-DD 
EXPSTART=        53105.9428153 /Start time of obs. sequence (MJD) 
EXPEND  =        53105.9637181 /End time of obs. sequence (MJD) 
MINWAVE =        1675.66137068 /Minimum wavelength 
MAXWAVE =        10198.7137562 /Maximum wavelength 
MAXFLUX =          9.14272E-13 /Maximum Flux in erg/sec/cm^2/Angstrom 
OFFSETPX=            -0.324710 /Offset of star from center of slit  
(pixels) 
 
 
OFFSETPX gives the offset of the target in pixels (0.05 arcsecond/pixel) from the center 
of the slit as measured using the G750L fringe flats (Section 3.2).  DATAQUAL is used 
to flag data that may be suspect because of a large centering error within the slit.  The 
data are flagged as suspect if the offset error is more than 0.9 pixels.  As center of the 
target star moves off the slit, the estimated offset computed from the fringe flat does not 
continue to increase and may even get smaller because the fringe flat is being correlated 
with the wings of the stellar PSF.  To identify this problem, we also flag data as suspect if 
the star’s visual magnitude determined from the observed spectrum differs by more than 
0.1 from the visual magnitude given in the Tycho II catalog or if the B-V magnitude 
differs by more than 0.3. 
 
FITQUAL gives an indicator of the quality of the fit of the stellar model to the spectrum: 
 1 – indicates that a good fit was found. 
 2 – indicates that we were able to f it the spectrum but the RMS (see section 4)  
was greater 0.035 
 3 – indicates that no fit was possi ble or the RMS was greater than 0.08 
 
If a fit was performed, keywords FIT_RMS, TEFF, LOG_G, LOG_L, ALPHA, 
EBMV, and DPC will be populated.  DPC is the estimated distance in parsecs based on 
the best fit model that falls on a Victoria-Regina isochrone with the restriction that the 
parallax for the distance falls within a 2-sigma error of new Hipparcos catalog parallax. 
 
The file’s binary table extension contains five columns: 
 
WAVELENGTH – wavelength in Angstroms 
DATAQUAL= 'good    '           /Data Quality: good or suspect 
HIPP    =                80822 /Hipparcos number 
FITQUAL =                    1 /Fit Quality: 1=good, 2=poor, 3=not 
fit 
FIT_RMS =            0.0220523 /RMS of FIT 
TEFF    =              5875.56 /Effective Temperature (K) 
LOG_G   =              2.77228 /Log Surface Gravity 
LOG_Z   =            -0.600000 /Log Abundance (with respect to Solar) 
ALPHA   = 'n       '           /"a" if alpha-enhanced 
EBMV    =             0.125778 / E(B-V) derived from fit 
DPC     =              606.252 /Distance (Parsecs) 
HISTORY STIS Next Generation Spectral Library, Version 1, 27-Feb-2008 
HISTORY  -------------------------------------------- 
HISTORY        FILENAME        OPT_ELEM   TEXPTIME 
HISTORY  -------------------------------------------- 
HISTORY  o8ru1p010_raw.fits   G230LB         600.000 
HISTORY  o8ru1p020_raw.fits   G230LB         600.000 
HISTORY  o8ru1pjiq_raw.fits   G430L          30.0000 
HISTORY  o8ru1pjjq_raw.fits   G430L          30.0000 
HISTORY  o8ru1pjkq_raw.fits   G750L          30.0000 
HISTORY  o8ru1pjlq_raw.fits   G750L          30.0000 
END 
 
 
FLUX – observed flux in ergs/cm2/second/Angstrom 
STATERR – propagated counting statistical errors. 
FLUX_UNRED – flux after extinction correction using the E(B-V) computed 
during model fitting (section 4) 
FLUX_10PC – Flux scaled to a stellar distance of 10 parsecs (absolute flux).  The 
distance was computed during the model fitting (section 4) allowing up to a 
two-sigma error in the parallax from the Hipparcos catalog. 
 
2. Improvement s to the Standard Pipeline Calibration 
 
The following improvements to the standard STIS pipeline processing were used to 
process the individual NGSL observations: 
 
1) New spectral trace files were produced used. These were constructed from the 
average spectral y-position versus x-position for the 52X0.2E1 aperture for all of 
the NGSL stars.   
2) Custom fringe flats were constructed for each of the G750L observations using 
the tungsten lamp observation associated with the visit. 
3) A more appropriate background subtraction was used.  The standard pipeline uses 
only a lower background taken 300 pixels below the spectrum for the E1 aperture 
position.  We used both an upper and lower background taken much closer to the 
spectrum (30 pixels). 
4) A larger extraction slit height (11 instead of 7) was used to improve overall 
photometric precision at the cost of some loss in S/N.  
5) Custom wavelength dispersion coefficients were created for the E1 aperture 
position. 
6) Wavecals were not taken with the stellar observations.  Zero point offsets were 
computed for each spectrum using the positions and wavelengths of strong stellar 
features. 
7) Custom sensitivity curves applicable to our background subtraction method were 
created for the 52X0.2E1 aperture using observations of BD+75D325 centered in 
the aperture. 
 
3. Post Pipeline Processing. 
 
After routine pipeline processing the following steps were used to construct the final 
composite spectra: 
1. The observations for each grating were obtained at two slightly different y-
positions along the slit.  These were co-added after rejection of cosmic rays and 
hot pixels.   
2. The observations of the three gratings were merged with an average taken in the 
overlap region (2990 – 3060 Angstroms for G230LB and G430L, 5500-5650 
Angstroms for G430L and G750L).   
3. In-order grating scatter was subtracted for G230LB (Section 3.1). 
4. A correction was performed for the slit throughput for targets not centered in the 
slit (Section 3.2) 
 
3.1  G230LB Scattered Light Correction 
 
There is significant in-order grating scatter which affects the flux values for G230LB.   
This scattered light is particularly evident for very cool stars.  We have modeled this 
scatter by: 
 
 SC = C0(1+Slope*x) 
 
Where SC is the scattered light in counts/second/pixel in the extracted net spectrum, X is 
the pixel number from 0 to 1023, C0 is the amount of scattered light, and Slope is the 
slope of the scattered light profile across the detector. 
 
To construct the model, we selected 103 red stars in the NGSL library with a B-V>1.1, 
Vmag<10, and the position of the target within 0.9 pixels of the center of the 52X0.2E1 
aperture. (as measured using the G750L fringe flat, see section 3.2).   
 
For each star, the scattered light was fit by a straight line between pixels 10 and 220 
where actual stellar flux should be negligible for the 103 red stars. 
 
 Scattered light  = c0 + c1*x       (1) 
 
 
Figure 1 (G230LB Slope of the scattered light) 
 
Figure 1 shows the values of c1 plotted versus c0.This was fit by a straight line with a 
slope of 0.0024. 
 
Thus equation 1 can be written as: 
 
 Scattered Light = c0 (1+0.0024*x)      (2) 
 
We model the amount of scatter, c0, by: 
 
  
n
CAc ∫= λλ /)(0  
 
Where C(λ) is the scattered light in counts/second and n is a selectable exponent.  The 
integration is done from 2000 to 10,000 Angstroms.  The limits of integration present a 
problem since we only get the net count rates below 3060 Angstroms for a G230LB 
observation.  To get the net count rates for G230LB beyond the portion visible on the 
detector we use the flux measured by the G430L and G750L gratings.  These can be 
converted to G230LB count rates using the combined G230L component efficiencies 
measured prior to launch.  Figure 2 shows the G230L sensitivity determined from the 
combined component efficiencies (solid line) compared to the post-launch G230LB 
sensitivity shown as a dashed line. 
 
Figure 2 (G230LB Sensitivity) 
 
We have chosen the exponent, n=3, to minimize the scatter in the fit for the 103 NGSL 
red stars.  Figure 3 shows a plot of C0 versus the integral of the count rate divided by 
wavelength cubed.  The slope of the fitted line gives a value of A=5523. 
 
 
Figure 3 (G230LB scattered light) 
 
 
Thus the final scattered light correction algorithm is given by: 
 
1. Process the observations for all three gratings using no scattered light correction. 
2. Merge the flux from the three gratings and divide by the G230LB sensitivity from 
2000 to 10000 Angstroms to get G230LB count rates. 
3. Integrate the count rate divided by the wavelength cubed from 2000 to 10000 
Angstroms and multiply by 5523 to get C0. 
4. Subtract C0*(1+0.0024)*X from the G230L net spectrum and multiply by the 
post-launch sensitivity curve. 
 
One concern is that the straight line fit to the scattered light which is very good for the 
first 220 pixels may not an accurate extrapolation to pixel 1023.  To test the extrapolation 
we can use the overlap region between the G230LB and G430L gratings.  Figure 4 shows 
the ratio of the G430L/G230LB  measured flux from 2980 to 3040 Angstroms versus B-
V when no G230LB scattered light correction is performed.  Stars with a B-V greater 
than 1 show the ratio becoming smaller as B-V increases indicating that the G230LB flux 
values are too large.  The same plot (Figure 5) when the scattered light correction is done 
shows much better results for these redder stars. 
 
  
 
 
Figure 4 (No G230LB scattered light correction) 
 
 
Figure 5 (With G230LB Scattered light correction) 
 
 
Figure 6 shows an example of the G230LB spectrum before (black) and after (green) 
scattered light correction.  For comparison, the red plot shows a STIS Echelle observation 
(E230M) binned to the same resolution. 
 
Figure 6: G230M scatter light correction compared to an E230M spectrum 
 
 
3.2 Correction for position of the target within the Slit. 
 
Significant photometric errors occur when the target star is not properly centered in the 
narrow 52X0.2 slit (as indicated in Figure 1).   Fortunately, we can measure the offset of 
the targets from the center of the slit using the G750L fringe flats taken with the 
observation.  If a target is centered in the slit, the fringes in the flat field observation will 
align with the fringes in the stellar observation.   If the target is not centered, the offset 
between the fringes in the flat and the stellar spectra give you the offset within the slit.  
 
Our first adjustment is for the change in the visual magnitude measured from the 
spectrum.   Figure 7 shows the difference in the VT magnitude computed from the NGSL 
spectra and those tabulated in the Tycho II catalog plotted versus target offset from the 
center of the slit.   A least square quadratic is fit to the difference (solid line in Figure 7) 
allows us to compute a scale factor to correct the observed spectral flux. 
 
Figure 7: NGSL magnitude errors resulting when target is not centered within the 
slit. 
 
Since the slit correction is wavelength dependent we investigated the B-V computed from 
the spectra as compared to the Tycho II values.  Figure 8 shows the small systematic 
differences as a function of the target’s position within the slit.  We fit the residuals with 
a least squares cubic spline curve (solid curve in Figure 8) and used the fit to correct the 
spectra.  The spectral flux is corrected with a correction factor (a linear function versus 
wavelength between 2900 and 5700 Angstroms  normalized to 1.0 at 5300 Angstroms.).  
The slope of the function is computed to adjust the B-V to the center of the slit using the 
solid correction curve in Figure 8.  The correction factor at 5700 Angstroms was used to 
correct all data above 5700 Angstroms. The 2900 Angstrom  correction factor was used 
for all data below 2900 Angstroms. 
 
 
Figure 8: NGSL Spectroscopic BT-VT compared to Tycho II catalog values. 
 
 
To correct wavelengths above 5700, we used best fit stellar models (using only the data 
below 5700 to determine the model parameters).  We compared the ratio of the calibrated 
NGSL spectra to the model for wavelengths above 5700 Angstroms.  Figure 9 shows the 
median ratio of the NGSL/Model spectra for various ranges of slit offsets.   Figure 10 
shows the flux ratios versus aperture position for all observations (plotted as diamonds) at 
various wavelengths (in 400 Angstrom bins).  We fit a smooth spline curve versus slit 
offset to determine a correction for each wavelength bin.  These corrections were then 
applied to the calibrated spectrum.  Figure 11 shows the same plot as Figure 10 after the 
data has been corrected 
 
 
Figure 9: Flux errors when target is not centered within the target aperture. 
 
 
Figure 10: Flux correction versus position within the aperture in 400 Angstrom 
wavelength bins 

 
 
Figure 11: Observed NGSL spectra compared to the Munari models after 
correction for the offset of the target within the aperture. 
 
As a final check of the slit correction we compared the ratio of the NGSL spectra with the 
matching spectra in the Miles ground base spectral library.   Figure 12 shows the ratio (in 
50 Angstrom bins) for various aperture offsets before correcting the NGSL spectra for the 
position of the star within the slit.  Figure 13 shows the same comparison after slit 
correction.  Results indicate that the slit correction is good to within a couple of percent 
and the overall sensitivity calibration of the Miles observations agree with the NGSL 
calibration to within approximately 3 percent from 3600 to 7300 Angstroms. 
 
 
Figure 12: Comparison of the STIS NGSL spectra (without a slit correction) with 
the ground based Miles spectra for various offsets from the center of the slit. 
 
 
 
Figure 13: Comparison of the STIS NGSL spectra (after the slit correction) with the 
ground based Miles spectra for various offsets from the center of the slit. 
 
4.0 Stellar Mod el Fitting 
 
We estimated the effective temperature, surface gravity, metallicity, and reddening of 
NGSL stars by comparing the observed spectrum to spectral models by Castelli (2004). 
The least-squares fit attempts to minimize the RMS difference between the model and 
observed spectrum (both normalized to an average of 1.0 in the range of 4500 to 7000 
Angstroms). The models were restricted to those falling on an Victoria-Regina isochrone 
(Vandenberg et al. 2006)  The isochrone surface gravity was derived on the assumption 
that the target distance is within a two-sigma error of the distance in the new Hipparcos 
catalog.  The following table gives the result for the best models.  The fit quality 
(FITQUAL keyword in the header) is set to 1(good) if the RMS is less than 0.035.  The 
quality is set to 2 (poor) if the RMS is between 0.035 and 0.08.  If the RMS is above 0.08 
or if no valid model good be found for the target (with the distance range from 
Hipparcos), no fit parameters are given.  
 
 
 Target        Data     Fit   FIT    Teff   Log(g) Log(Z) alpha 
  Name       Quality Quality  RMS 
BD+112998      good      1   0.022   5876.    2.8   -0.6    n 
BD+363168      suspect   3 
BD+413306      good      1   0.027   5156.    4.5   -0.3    n 
BD-122669      good      1   0.018   6960.    4.1   -1.9    a 
BD092860       suspect   1   0.030   5926.    3.5   -0.7    a 
BD174708       good      1   0.018   6174.    3.8   -1.9    a 
BD292091       suspect   1   0.017   5803.    4.6   -2.0    a 
BD371458       suspect   1   0.015   5463.    3.1   -2.0    a 
BD413931       suspect   1   0.019   5486.    4.7   -1.8    a 
BD423607       suspect   1   0.022   5963.    4.5   -2.0    a 
BD442051       good      2   0.072   3741.    4.9   -0.3    n 
BD511696       good      1   0.020   5776.    4.5   -1.3    a 
BD592723       good      1   0.020   6419.    4.4   -2.0    a 
BD660268       suspect   1   0.030   5453.    4.7   -1.8    a 
BD720094       good      1   0.016   6286.    3.9   -2.0    a 
CD-259286      suspect   1   0.018   6426.    3.0   -1.2    a 
CD-3018140     suspect   1   0.015   6334.    3.9   -2.0    a 
CD-621346      good      1   0.023   5652.    3.2   -0.6    n 
CD-691618      good      3 
G019-013       good      1   0.027   4452.    4.6   -0.3    n 
G021-024       suspect   2   0.037   4230.    4.7   -0.1    n 
G029-023       good      1   0.019   6176.    3.7   -2.0    a 
G114-26        good      1   0.019   6144.    4.1   -1.8    a 
G115-58        suspect   1   0.023   6156.    3.7   -1.6    a 
G12-21         suspect   1   0.020   6034.    4.1   -1.4    a 
G13-35         good      1   0.016   6164.    3.9   -1.8    a 
G169-28        suspect   1   0.019   5749.    3.9   -1.5    a 
G17-25         suspect   1   0.022   5289.    4.6   -1.0    a 
G18-39         good      1   0.017   5952.    3.6   -1.7    a 
G18-54         good      1   0.022   5920.    3.8   -1.7    a 
G180-24        good      1   0.018   6109.    4.0   -1.5    a 
G187-40        suspect   1   0.017   5776.    4.0   -1.6    a 
G188-22        suspect   1   0.016   5926.    3.9   -1.5    a 
G188-30        suspect   1   0.030   5497.    4.6   -1.2    a 
G192-43        suspect   1   0.017   6163.    3.9   -1.6    a 
G194-22        suspect   1   0.018   6193.    4.1   -1.6    a 
G196-48        suspect   1   0.013   5704.    3.4   -1.9    a 
G20-15         suspect   1   0.017   6089.    4.1   -1.8    a 
G202-65        suspect   1   0.028   6989.    4.4   -1.0    a 
G231-52        suspect   1   0.023   5426.    4.7   -2.0    a 
G234-28        good      1   0.019   6097.    3.9   -1.6    a 
G24-3          good      1   0.018   6053.    3.9   -1.8    a 
G243-62        suspect   1   0.028   4693.    4.8   -1.7    a 
G260-36        suspect   1   0.026   5057.    4.5    0.1    n 
G262-14        suspect   1   0.025   5085.    3.8   -0.7    a 
G63-26         suspect   1   0.024   6123.    3.9   -1.6    a 
G88-27         good      1   0.018   6141.    3.7   -1.6    a 
GJ825          good      2   0.054   4015.    4.7    0.1    n 
GL109          suspect   3 
GL15B          good      2   0.079   3741.    4.9   -0.2    n 
HD000319       good      1   0.021   8195.    3.9   -0.4    n 
HD000886       good      3 
HD001461       good      1   0.026   5790.    4.3    0.3    n 
HD002665       good      1   0.011   5082.    2.2   -2.0    a 
HD002857       suspect   2   0.036   7708.    2.8   -0.6    n 
HD003360       good      3 
HD003712       suspect   1   0.022   4603.    1.8   -0.3    n 
HD004128       suspect   1   0.021   4815.    2.3   -0.2    n 
HD004727       good      3 
HD004813       good      1   0.024   6252.    4.3   -0.1    n 
HD005256       good      1   0.021   5176.    3.6   -0.6    a 
HD005395       good      1   0.018   4808.    2.5   -0.5    n 
HD005544       good      1   0.016   4474.    1.7   -0.5    n 
HD005916       good      1   0.017   4926.    2.5   -0.6    n 
HD006229       good      1   0.017   5489.    2.5   -0.6    n 
HD006734       good      1   0.021   4904.    3.1   -0.6    a 
HD006755       good      1   0.013   5064.    2.6   -2.0    a 
HD008491       good      1   0.019   4688.    2.5   -0.1    n 
HD008724       suspect   1   0.019   4552.    1.3   -2.0    a 
HD008890       good      3 
HD009051       suspect   1   0.012   4842.    2.0   -1.8    a 
HD010380       good      1   0.015   4063.    1.2   -0.5    n 
HD010780       good      1   0.024   5357.    4.4    0.1    n 
HD012533       suspect   1   0.020   4225.    1.3   -0.2    n 
HD013520       suspect   1   0.018   3920.    0.8   -0.6    n 
HD015089       suspect   2   0.039   8823.    4.2    0.5    n 
HD016031       suspect   1   0.017   6299.    3.9   -1.8    a 
HD017072       good      1   0.020   5709.    3.1   -0.6    a 
HD017081       good      1   0.029  13057.    3.6   -0.5    n 
HD017361       good      1   0.019   4603.    2.4   -0.1    n 
HD017925       good      1   0.022   5208.    4.5    0.1    n 
HD018078       good      2   0.037   7404.    2.9    0.5    n 
HD018769       good      1   0.028   8382.    4.1    0.5    n 
HD018907       suspect   1   0.025   5309.    3.7   -0.1    n 
HD019019       good      1   0.023   6197.    4.4    0.0    n 
HD019308       good      1   0.022   5842.    4.1    0.2    n 
HD019445       good      1   0.018   6200.    4.3   -1.9    a 
HD019656       suspect   1   0.017   4563.    2.1   -0.3    n 
HD019787       good      1   0.020   4771.    2.6   -0.1    n 
HD020039       good      1   0.019   5200.    3.5   -0.6    a 
HD020630       good      1   0.026   5687.    4.4    0.0    n 
HD021742       good      1   0.024   5342.    4.4    0.4    n 
HD022049       good      1   0.021   5130.    4.5    0.0    n 
HD022484       good      1   0.023   6141.    4.1    0.1    n 
HD023439       good      1   0.021   5249.    4.5   -0.6    a 
HD025329       good      1   0.028   4930.    4.7   -1.1    a 
HD025893       good      1   0.023   5342.    4.4    0.1    n 
HD025975       good      1   0.020   4942.    3.3    0.0    n 
HD026297       good      1   0.012   4423.    1.0   -2.0    a 
HD026630       good      2   0.056   4542.    1.1   -2.0    a 
HD027295       good      2   0.058  11548.    4.2   -0.4    n 
HD028946       good      1   0.024   5389.    4.4   -0.0    n 
HD028978       suspect   1   0.023   8749.    3.5    0.1    n 
HD029391       good      1   0.027   7522.    4.3    0.2    n 
HD029574       good      1   0.025   4049.    0.6   -1.7    a 
HD030614       good      3 
HD030834       good      1   0.017   4008.    0.9   -0.6    n 
HD031219       suspect   1   0.023   6087.    3.9    0.3    n 
HD031421       good      1   0.018   4434.    1.9   -0.5    n 
HD033793       good      3 
HD034078       good      3 
HD034797       suspect   2   0.036  12884.    4.1    0.0    n 
HD034816       good      3 
HD036702       good      1   0.014   4319.    0.8   -2.0    a 
HD036960       suspect   3 
HD037202       good      3 
HD037216       suspect   1   0.023   5437.    4.4   -0.0    n 
HD037763       good      1   0.020   4597.    2.9    0.2    n 
HD037828       good      1   0.012   4359.    1.0   -1.6    a 
HD038237       good      1   0.025   8100.    3.9    0.1    n 
HD038510       good      1   0.019   5897.    4.0   -0.8    a 
HD039587       good      1   0.026   6057.    4.4    0.1    n 
HD039833       good      1   0.027   5926.    4.3    0.2    n 
HD040573       good      1   0.018  10200.    4.2   -0.4    n 
HD041357       good      1   0.025   7578.    3.5    0.2    n 
HD041661       good      1   0.019   6360.    3.2   -0.1    n 
HD041667       suspect   1   0.019   4623.    1.8   -1.0    a 
HD043042       suspect   1   0.022   6542.    4.2    0.0    n 
HD044007       good      1   0.014   4930.    2.3   -1.8    a 
HD045282       good      1   0.015   5199.    2.8   -1.8    a 
HD046703       good      2   0.046   6484.    2.5   -0.0    n 
HD047839       suspect   3 
HD048279       good      3 
HD050420       good      1   0.022   7167.    2.8   -0.0    n 
HD052973       suspect   3 
HD055057       good      1   0.022   7289.    3.5    0.2    n 
HD055496       suspect   1   0.024   4978.    2.1   -0.6    n 
HD057060       suspect   3 
HD057061       suspect   3 
HD057727       good      1   0.020   4942.    2.9   -0.3    n 
HD058343       suspect   3 
HD058551       good      1   0.020   6311.    4.2   -0.4    n 
HD059612       good      3 
HD060319       good      1   0.018   5926.    4.0   -0.9    a 
HD061064       good      1   0.021   6600.    3.4    0.2    n 
HD061603       good      1   0.027   3841.    0.8   -0.3    n 
HD062412       good      1   0.021   4876.    2.6   -0.1    n 
HD063077       good      1   0.021   5926.    4.2   -0.7    a 
HD063700       good      1   0.027   4498.    1.5   -0.5    n 
HD063791       good      1   0.011   4778.    1.6   -2.0    a 
HD064412       good      1   0.020   5742.    4.1   -0.6    a 
HD065228       suspect   1   0.021   5926.    2.4    0.4    n 
HD065354       good      1   0.021   3897.    0.8   -0.6    n 
HD065714       good      1   0.023   4882.    2.2    0.0    n 
HD067390       suspect   1   0.022   7193.    3.6   -0.2    n 
HD068988       good      1   0.026   5897.    4.3    0.3    n 
HD071160       good      1   0.017   3949.    0.9   -0.5    n 
HD072184       good      1   0.019   4620.    2.6    0.1    n 
HD072324       good      1   0.020   4793.    2.3   -0.1    n 
HD072505       good      1   0.019   4519.    2.3    0.1    n 
HD072968       suspect   2   0.066   9253.    3.9    0.5    n 
HD073710       good      1   0.021   4897.    2.2    0.1    n 
HD074088       good      1   0.028   3740.    0.5   -1.1    a 
HD074721       good      1   0.027   8774.    3.3   -0.6    n 
HD076291       good      1   0.017   4534.    2.3   -0.2    n 
HD076932       good      1   0.019   6034.    4.1   -0.7    a 
HD078316       good      2   0.073  12442.    3.7   -0.1    n 
HD078362       good      1   0.028   7120.    3.8    0.5    n 
HD078479       good      1   0.018   4476.    2.1    0.1    n 
HD079158       good      3 
HD079349       suspect   1   0.027   3834.    1.1   -0.1    n 
HD079469       good      1   0.027  10489.    4.2   -0.2    n 
HD080607       good      1   0.023   5511.    3.7    0.4    n 
HD081797       suspect   1   0.018   4188.    1.4   -0.1    n 
HD082395       good      1   0.018   4664.    2.4   -0.3    n 
HD082734       good      1   0.022   4949.    2.4    0.2    n 
HD083212       good      1   0.014   4448.    1.2   -1.8    a 
HD085380       good      1   0.024   6142.    4.0    0.2    n 
HD086322       good      1   0.019   4687.    2.3   -0.3    n 
HD086986       good      1   0.029   8230.    3.5   -1.2    a 
HD087140       good      1   0.011   5093.    2.3   -2.0    a 
HD087737       good      3 
HD090862       good      1   0.017   3937.    0.9   -0.6    n 
HD091316       good      3 
HD093329       good      1   0.021   8323.    3.3   -0.6    n 
HD093813       good      1   0.015   4264.    1.6   -0.5    n 
HD094028       good      1   0.017   6104.    4.2   -1.4    a 
HD095241       good      1   0.021   6034.    3.7   -0.1    n 
HD095735       good      3 
HD095849       good      1   0.018   4459.    1.9   -0.0    n 
HD096446       good      3 
HD097633       good      1   0.020   9107.    3.6   -0.2    n 
HD099648       good      1   0.019   4811.    2.0   -0.3    n 
HD101013       good      1   0.023   4837.    2.5   -0.0    n 
HD101107       good      1   0.022   6937.    3.9   -0.1    n 
HD102212       good      2   0.042   3800.    1.1    0.1    n 
HD102780       good      1   0.032   3669.    0.5   -0.6    n 
HD103036       suspect   1   0.013   4282.    0.8   -1.9    a 
HD105452       good      1   0.023   7120.    4.2   -0.2    n 
HD105546       good      1   0.016   5123.    2.3   -1.8    a 
HD105740       good      1   0.016   4642.    2.2   -0.7    a 
HD106304       suspect   1   0.015   9376.    3.6   -1.8    a 
HD106516       good      1   0.020   6310.    4.3   -0.7    a 
HD107582       good      1   0.022   5642.    4.3   -0.6    a 
HD108945       suspect   1   0.033   8642.    3.9    0.5    n 
HD109387       good      3 
HD109995       good      1   0.015   8251.    3.3   -1.6    a 
HD110073       good      2   0.069  12041.    3.8   -0.4    n 
HD110885       good      1   0.017   5823.    2.9   -0.8    a 
HD111464       good      1   0.017   4004.    1.0   -0.6    n 
HD111515       good      1   0.025   5549.    4.4   -0.3    n 
HD111721       good      1   0.013   4911.    2.2   -1.9    a 
HD111786       good      1   0.022   7598.    3.9   -1.3    a 
HD112413       good      2   0.051  11687.    4.0    0.5    n 
HD113002       good      1   0.020   5342.    2.7   -0.6    n 
HD113092       good      1   0.015   4203.    1.3   -0.6    n 
HD114330       suspect   1   0.021   9289.    3.6   -0.1    n 
HD114710       good      1   0.026   6142.    4.4    0.2    n 
HD115617       good      1   0.026   5557.    4.3    0.0    n 
HD117880       suspect   1   0.029   9426.    3.7   -0.6    n 
HD118055       good      1   0.014   4320.    0.9   -1.8    a 
HD119971       good      1   0.019   4034.    1.1   -0.8    a 
HD121146       good      1   0.017   4387.    2.2   -0.1    n 
HD122064       good      1   0.021   4908.    4.5    0.4    n 
HD122956       good      1   0.011   4625.    1.4   -2.0    a 
HD123657       good      3 
HD124186       good      1   0.019   4420.    2.2    0.2    n 
HD124425       good      1   0.020   6389.    3.6   -0.0    n 
HD124547       suspect   1   0.016   4084.    1.2   -0.4    n 
HD126327       suspect   3 
HD126511       good      1   0.025   5486.    4.3    0.2    n 
HD126614       good      1   0.027   5437.    4.0    0.5    n 
HD126661       good      1   0.028   7757.    3.6    0.4    n 
HD128000       good      1   0.025   3837.    0.9   -0.4    n 
HD128279       good      1   0.013   5456.    3.0   -2.0    a 
HD128801       good      1   0.014  10123.    3.7   -1.9    a 
HD128987       suspect   1   0.025   5537.    4.4   -0.0    n 
HD131873       good      1   0.017   3926.    0.9   -0.6    n 
HD132345       good      1   0.023   4442.    2.3    0.2    n 
HD132475       good      1   0.016   5687.    3.6   -1.6    a 
HD134113       good      1   0.022   5820.    3.9   -0.6    a 
HD134439       suspect   1   0.032   5130.    4.7   -1.5    a 
HD134440       good      1   0.028   4970.    4.7   -1.1    a 
HD136726       good      1   0.016   4093.    1.3   -0.4    n 
HD137759       suspect   1   0.028   4563.    2.3   -0.2    n 
HD137909       suspect   2   0.045   7663.    3.8    0.0    n 
HD138716       suspect   1   0.019   4771.    2.9   -0.1    n 
HD138749       suspect   2   0.066  14141.    3.8   -0.6    n 
HD140232       good      1   0.028   8157.    4.2    0.5    n 
HD141795       good      1   0.026   8419.    4.2    0.3    n 
HD141851       good      1   0.021   8231.    4.0   -0.2    n 
HD142091       good      1   0.020   4749.    2.9   -0.1    n 
HD142703       good      1   0.017   7337.    3.9   -1.4    a 
HD142860       suspect   1   0.023   6397.    4.2   -0.1    n 
HD142926       good      1   0.035  11908.    3.8    0.1    n 
HD143107       suspect   1   0.015   4303.    1.7   -0.4    n 
HD143459       good      1   0.023   9878.    3.6   -0.6    n 
HD145328       suspect   1   0.020   4764.    2.8   -0.2    n 
HD146051       good      2   0.040   3667.    0.7   -0.3    n 
HD146233       suspect   1   0.024   5842.    4.4    0.1    n 
HD147394       good      3 
HD147550       good      1   0.033  10074.    3.9   -0.0    n 
HD148293       good      1   0.021   4642.    2.4    0.0    n 
HD148513       good      1   0.021   4003.    1.3   -0.3    n 
HD149161       good      1   0.022   3884.    1.1   -0.3    n 
HD149382       good      3 
HD155763       good      3 
HD156283       good      1   0.018   4133.    1.3   -0.1    n 
HD157244       good      1   0.027   4206.    1.3   -0.1    n 
HD159181       good      2   0.080   4542.    1.5   -0.6    n 
HD160346       good      1   0.023   5004.    4.5    0.1    n 
HD160762       good      3 
HD160922       good      1   0.023   6600.    4.0    0.0    n 
HD161770       good      1   0.017   5818.    3.8   -1.6    a 
HD163346       good      2   0.074   5597.    2.5   -0.6    n 
HD163641       good      2   0.041  11482.    3.9   -0.1    n 
HD163810       suspect   1   0.019   5876.    4.4   -1.2    a 
HD164058       good      1   0.024   3896.    0.9   -0.4    n 
HD164257       good      2   0.043   7977.    3.5   -0.1    n 
HD164353       good      3 
HD164402       good      3 
HD164967       good      1   0.026   8534.    4.1   -0.6    n 
HD165195       good      1   0.034   4060.    0.6   -1.8    a 
HD165341       good      1   0.022   5319.    4.4    0.1    n 
HD166208       good      1   0.025   5107.    2.3    0.0    n 
HD166229       good      1   0.020   4576.    2.5    0.1    n 
HD166283       good      1   0.028   8289.    4.1    0.1    n 
HD166991       good      1   0.021   8497.    4.0   -0.3    n 
HD167006       good      2   0.059   3536.    0.4   -0.3    n 
HD167105       suspect   1   0.031   8574.    3.4   -0.6    n 
HD167278       suspect   1   0.023   6609.    4.2   -0.1    n 
HD167946       good      1   0.030  10634.    4.3   -0.1    n 
HD169191       good      1   0.015   4245.    1.6   -0.4    n 
HD170737       good      1   0.017   4942.    3.0   -0.9    a 
HD170756       good      2   0.040   6252.    2.4   -0.0    n 
HD170973       good      2   0.064  10676.    3.3    0.5    n 
HD172230       good      1   0.030   7542.    3.4    0.5    n 
HD172506       suspect   1   0.027   7226.    4.1    0.0    n 
HD173158       good      2   0.038   3856.    0.5   -1.3    a 
HD173819       suspect   2   0.071   3714.    0.6   -0.6    a 
HD174240       good      1   0.019   9274.    3.8   -0.2    n 
HD174959       suspect   1   0.030  14123.    3.7   -0.6    n 
HD174966       suspect   1   0.025   7693.    4.0    0.1    n 
HD175156       good      3 
HD175305       suspect   1   0.016   5037.    2.5   -1.6    a 
HD175545       good      1   0.018   4420.    2.4    0.0    n 
HD175640       good      2   0.057  11237.    3.9   -0.2    n 
HD175674       good      1   0.022   4242.    1.6   -0.1    n 
HD175805       good      1   0.021   6274.    3.6    0.2    n 
HD175865       suspect   3 
HD176232       suspect   1   0.030   7615.    3.9    0.3    n 
HD176437       good      3 
HD181720       good      1   0.021   5709.    3.9   -0.6    a 
HD183324       suspect   1   0.024   8673.    4.1   -1.1    a 
HD183915       good      1   0.022   4319.    1.4   -0.4    n 
HD184266       suspect   1   0.018   6109.    3.0   -1.1    a 
HD185144       good      1   0.022   5293.    4.5   -0.1    n 
HD185351       good      1   0.020   4904.    3.0   -0.1    n 
HD187111       good      1   0.015   4423.    1.1   -1.8    a 
HD187879       good      3 
HD188262       suspect   3 
HD190073       good      3 
HD190360       good      1   0.027   5531.    4.2    0.2    n 
HD190404       good      1   0.026   5110.    4.6   -0.3    n 
HD191026       good      1   0.022   5149.    3.7    0.1    n 
HD191277       good      1   0.018   4476.    2.5    0.1    n 
HD193281       good      1   0.018   8249.    3.6   -0.6    n 
HD193495       good      3 
HD194093       suspect   3 
HD194453       suspect   1   0.023  10241.    3.9    0.0    n 
HD195434       suspect   1   0.026   5152.    4.6   -0.5    n 
HD196218       good      1   0.023   6367.    4.2    0.0    n 
HD196426       good      1   0.029  12800.    4.0   -0.5    n 
HD196662       suspect   1   0.031  14204.    3.7   -0.6    n 
HD196725       good      1   0.022   4003.    0.9   -0.6    n 
HD196892       good      1   0.018   6086.    4.1   -1.0    a 
HD197177       good      1   0.019   4748.    1.9   -0.4    n 
HD198809       good      1   0.020   5291.    2.9   -0.0    n 
HD200081       good      1   0.034   5120.    2.6   -0.5    n 
HD200905       good      1   0.030   3793.    0.6   -0.6    n 
HD201091       good      1   0.030   4404.    4.6   -0.1    n 
HD201377       good      1   0.023   8084.    3.9    0.0    n 
HD201601       good      1   0.034   7631.    4.0   -0.1    n 
HD203638       good      1   0.020   4520.    2.3   -0.1    n 
HD204041       good      1   0.020   8260.    4.2   -0.6    n 
HD204155       good      1   0.019   5787.    3.9   -0.6    a 
HD204543       suspect   1   0.014   4667.    1.4   -2.0    a 
HD204867       suspect   3 
HD205202       good      1   0.019   6564.    3.6   -0.6    a 
HD205811       good      1   0.021   9287.    4.2   -0.1    n 
HD206778       good      1   0.030   4095.    1.3    0.1    n 
HD210745       good      1   0.033   4003.    0.9   -0.6    n 
HD210807       good      1   0.018   4926.    2.3   -0.3    n 
HD212516       good      2   0.044   3578.    0.5   -0.5    n 
HD212593       good      3 
HD215665       good      1   0.019   4752.    2.0   -0.2    n 
HD217107       good      1   0.027   5631.    4.2    0.4    n 
HD217357       good      2   0.042   4153.    4.7   -0.3    n 
HD221377       good      1   0.020   6556.    3.9   -0.6    n 
HD222404       good      1   0.022   4790.    3.1    0.1    n 
HD224801       good      2   0.062  11757.    3.9    0.5    n 
HD224926       good      2   0.042  13499.    3.8   -0.5    n 
HD232078       suspect   2   0.072   3640.    0.5   -0.6    n 
HD284248       suspect   1   0.016   6171.    4.0   -1.8    a 
HD345957       good      1   0.017   5920.    3.8   -1.4    a 
HR0753         suspect   1   0.030   4882.    4.7   -0.2    n 
HR8086         suspect   2   0.047   4230.    4.7   -0.7    a 
VBNVUL         suspect   2   0.046   6937.    2.4   -0.5    n 
VGKCOM         suspect   3 
VIWCOM         good      3 
 
 