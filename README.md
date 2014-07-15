phot-sex
========

Matches stars in the field, constructs growth curve for stars in all bands, constructs SEDs for stars and galaxy, if given coordinates 
just uses catalogs generated from SExtractor

all results assume you use accompanying .sex and .sh files to generate data from images

Also builds a catalog of all sources in the field, converts all fluxes to uJy, and applies extinction corrections based on data from NED

**Columns of catalog (from SExtractor): (Note: all fluxes are ZP, extinction, and aperture corrected (assuming max flux at 2" and capping measured flux at 0.5"))**

NUMBER: Object # from SExtractor; we just match up rows between band catalogs, so this will break silently if the catalogs are not generated first with F160, then using F160 as a reference

X_IMAGE & Y_IMAGE: pixel coordinates of object, handy if you have the original images

X_WORLD & Y_WORLD: world coordinates of object, generated by SExtractor

FLUX_AUTO_160: flux in uJy of object, measured by SExtractor (corrected as above, minus aperture correction)

FLUXERR_AUTO_160: error in FLUX_AUTO_160

PHOT_RADIUS_1 & PHOT_RADIUS_2: fraction-of-light radii 0.5 & 0.9, measured in pixels

FLUX_APER_<BAND>: flux in uJy of the object (calculated with COGs and corrected as bove)

FLUXERR_APER_<BAND>: error in FLUX_APER_<BAND>, also in uJy
