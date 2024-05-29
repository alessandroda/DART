TROPOMI S5P data
========================

Overview
--------

The observations are the TROPOMI S5P satellite products (spatial resolution 5.5 x 3.5 km2 and daily coverage), more
specifically:

#. L2__NO2___, PUM-NO2 Nitrogen Dioxide (NO2), total and tropospheric columns
#. L2__SO2___, Sulfur Dioxide (SO2) from COBRA data processing (SO2CBR)
#. L2__HCHO___, Formaldehyde total column 

The file are those created by the Copernicus Satellite Operator code (CSO). `https://ci.tno.nl/gitlab/cams/cso`

Programs
--------

The programs in the ``DART/observations/obs_converters/S5P_TROPOMI_L3`` directory extract data from a file in 
convert_s5p_tropomi_l3.f90 and converts this into the observation sequence (obs_seq) files. Build them in the ``work`` 
directory by running the ``./quickbuild.sh`` script. In addition to the converters, the ``advance_time`` and ``obs_sequence_tool``
utilities will be built. 

There are currently converters for these data types:

=========================== ======================
S5P tropomi         convert_s5p_tropomi_l3
=========================== ======================

Example data files are in the ``data`` directory. Example scripts for converting batches of these files are in the
``shell_scripts`` directory. These are *NOT* intended to be turnkey scripts; they will certainly need to be customized
for your use. There are comments at the top of the scripts saying what options they include, and should be commented
enough to indicate where changes will be likely to need to be made.

The converter has a hard-coded input filename:

======================= ================= ================
convert_s5p_tropomi_l3: S5p_NO2_16685.nc  obs_seq.out   
======================= ================= ================
