import numpy as np
from netCDF4 import Dataset
from datetime import datetime

def convert_time_ref(year, month, day, hour, minute, second, ref_year):
    # Convert date and time to a reference time (in seconds since a reference year)
    ref_date = datetime(ref_year, 1, 1, 0, 0, 0)
    target_date = datetime(year, month, day, hour, minute, second)
    delta = target_date - ref_date
    return delta.total_seconds()

def tropomi_no2_total_col_extract(filein, fileout, file_pre, cwyr_mn, cwmn_mn, cwdy_mn, cwhh_mn, cwmm_mn, cwss_mn, cwyr_mx, cwmn_mx, cwdy_mx, cwhh_mx, cwmm_mx, cwss_mx, path_mdl, file_mdl, cnx_mdl, cny_mdl):
    wyr_mn = int(cwyr_mn)
    wmn_mn = int(cwmn_mn)
    wdy_mn = int(cwdy_mn)
    whh_mn = int(cwhh_mn)
    wmm_mn = int(cwmm_mn)
    wss_mn = int(cwss_mn)
    wyr_mx = int(cwyr_mx)
    wmn_mx = int(cwmn_mx)
    wdy_mx = int(cwdy_mx)
    whh_mx = int(cwhh_mx)
    wmm_mx = int(cwmm_mx)
    wss_mx = int(cwss_mx)
    nx_mdl = int(cnx_mdl)
    ny_mdl = int(cny_mdl)

    # Constants
    Ru = 8.316
    Rd = 286.9
    eps = 0.61
    molec_wt_no2 = 0.0480
    molec_wt_so2 = 0.0641
    AvogN = 6.02214e23
    msq2cmsq = 1.e4
    P_std = 1013.25
    grav = 9.8
    cone_fac = 0.715567
    du2molpm2 = 4.4615e-4
    du2molcpm2 = 2.6867e20

    day_secs_beg = whh_mn * 60. * 60. + cwmm_mn * 60. + cwss_mn
    day_secs_end = whh_mx * 60. * 60. + cwmm_mx * 60. + cwss_mx

    # Read model grid
    with Dataset(f"{path_mdl}/{file_mdl}", "r") as f:
        lon_mdl = f.variables["XLONG"][:]
        lat_mdl = f.variables["XLAT"][:]
        delx = f.getncattr("DX")
        cen_lat = f.getncattr("CEN_LAT")
        cen_lon = f.getncattr("CEN_LON")
        if cen_lon < 0:
            cen_lon += 360.
        truelat1 = f.getncattr("TRUELAT1")
        truelat2 = f.getncattr("TRUELAT2")
        moad_cen_lat = f.getncattr("MOAD_CEN_LAT")
        stand_lon = f.getncattr("STAND_LON")
        pole_lat = f.getncattr("POLE_LAT")
        pole_lon = f.getncattr("POLE_LON")

    # Process satellite data
    with open(fileout, "w") as fid:
        # Handle each file in the list
        file_list_a = !ls -1 {filein}*
        file_list_b = file_list_a[0].split()
        nfile = len(file_list_b)

        for ifile in range(nfile):
            file_in = file_list_b[ifile]
            if not file_in:
                continue
            
            indx = file_in.find(file_pre) - 1
            if indx == -1:
                continue

            file_str_yy = int(file_in[indx + 21: indx + 24])
            file_str_mm = int(file_in[indx + 25: indx + 26])
            file_str_dd = int(file_in[indx + 27: indx + 28])
            file_str_hh = int(file_in[indx + 30: indx + 31])
            file_str_mn = int(file_in[indx + 32: indx + 33])
            file_str_ss = int(file_in[indx + 34: indx + 35])
            file_end_yy = int(file_in[indx + 37: indx + 40])
            file_end_mm = int(file_in[indx + 41: indx + 42])
            file_end_dd = int(file_in[indx + 43: indx + 44])
            file_end_hh = int(file_in[indx + 46: indx + 47])
            file_end_mn = int(file_in[indx + 48: indx + 49])
            file_end_ss = int(file_in[indx + 50: indx + 51])
            file_str_secs = file_str_hh * 60. * 60. + file_str_mn * 60. + file_str_ss
            file_end_secs = file_end_hh * 60. * 60. + file_end_mn * 60. + file_end_ss
            
            if file_str_secs > day_secs_end or file_end_secs < day_secs_beg:
                continue
            
            print(f"{ifile} {file_in}")
            print(f"file str {file_str_secs} cycle end {day_secs_end}")
            print(f"file end {file_end_secs} cycle str {day_secs_beg}")

            # Read TROPOMI data
            with Dataset(file_in, "r") as tropomi_data:
                scanline = np.max(tropomi_data.variables["PRODUCT_scanline"])
                pixel = np.max(tropomi_data.variables["PRODUCT_ground_pixel"])

                for ilin in range(1, scanline + 1):
                    date_str = str(tropomi_data.variables["PRODUCT_time_utc"][ilin - 1])
                    yyyy_tropomi = int(date_str[0:4])
                    mn_tropomi = int(date_str[5:7])
                    dy_tropomi = int(date_str[8:10])
                    hh_tropomi = int(date_str[11:13])
                    mm_tropomi = int(date_str[14:16])
                    ss_tropomi = int(date_str[17:19])
                    tropomidate = convert_time_ref(yyyy_tropomi, mn_tropomi, dy_tropomi, hh_tropomi, mm_tropomi, ss_tropomi, 2010)
                    
                    if tropomidate < windate_min or tropomidate > windate_max:
                        continue

                    for ipxl in range(1, pixel + 1):
                        if tropomi_data
