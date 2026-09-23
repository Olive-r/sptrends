Original dataset
----------------
Dataset:
Blended Vegetation Health Product (Blended-VHP)

Provider:
NOAA Center for Satellite Applications and Research (STAR)

Reference:
NOAA Center for Satellite Applications and Research (STAR).
Blended Vegetation Health Product (Blended-VHP).
https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php

Original dataset characteristics
--------------------------------
Variable:
SMN (Smoothed Normalized Difference Vegetation Index)

Original temporal resolution:
Weekly

Original spatial resolution:
4 km

Original temporal coverage:
1982-present

Original spatial coverage:
Global land areas

Original variables available:
ND (raw NDVI and brightness temperature),
SM (smoothed NDVI and brightness temperature),
VH (Vegetation Condition Index, Temperature Condition Index
and Vegetation Health Index)

Processing
----------
The example dataset included with sptrends was derived from the original
Blended Vegetation Health Product (Blended-VHP) using the following
processing steps:

1. Annual mean NDVI calculated from weekly SMN observations.
2. Temporal subset extracted (1982-2023).
3. Reprojected to the Eckert IV equal-area projection.
4. Resampled to a regular 100 km × 100 km grid using bilinear interpolation.
5. Exported as LZW-compressed GeoTIFF files.

Purpose
-------
The purpose of this dataset is to provide a compact, globally consistent
example for demonstrating the functionality of the sptrends package.
Users requiring the original spatial resolution, weekly observations,
or the complete Vegetation Health Product should download the official
NOAA dataset.

Acknowledgement
---------------
Please acknowledge the NOAA Center for Satellite Applications and
Research (STAR) as the original source of the vegetation data.

Further information
-------------------
Always consult the official NOAA STAR documentation before using the
original dataset for scientific research.

NOAA STAR Vegetation Health Products:
https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php

Blended-VHP data repository:
https://www.star.nesdis.noaa.gov/data/pub0018/VHPdata4users/data/Blended_VH_4km/geo_TIFF/