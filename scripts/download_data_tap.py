import splusdata
from astroquery.gaia import Gaia
import warnings

import numpy as np
import time
import warnings

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import hstack, vstack, Table, join
from astropy.time import Time

def inner_merge_tables(
    t1, 
    t2, 
    t1_coord_cols=["ra", "dec"], 
    t2_coord_cols=["ALPHA_J2000", "DELTA_J2000"], 
    sep=1,
    add_metadata=True,
    match_prefix=""
):
    """
    Merge two tables based on sky coordinate matching, keeping only matched rows.
    
    This function performs an inner join between two tables based on spatial proximity
    of celestial coordinates. Only rows that have a match within the specified separation
    limit are included in the output.
    
    Parameters
    ----------
    t1 : astropy.table.Table
        First table (reference table)
    t2 : astropy.table.Table
        Second table (table to match)
    t1_coord_cols : list, optional
        Column names for RA and Dec in the first table, by default ["ra", "dec"]
    t2_coord_cols : list, optional
        Column names for RA and Dec in the second table, by default ["ALPHA_J2000", "DELTA_J2000"]
    sep : float, optional
        Maximum separation in arcseconds for a valid match, by default 1
    t1_epoch : str, optional
        Epoch of coordinates in first table, by default "J2000"
    add_metadata : bool, optional
        Whether to add merge metadata to the output table, by default True
    match_prefix : str, optional
        Prefix for match quality columns added to the output, by default "match_"
    
    Returns
    -------
    astropy.table.Table
        Merged table containing only rows that match between the two input tables
    
    Notes
    -----
    - The matching is performed using the astropy.coordinates.SkyCoord.match_to_catalog_sky method
    - Column name conflicts are handled by appending '_1' and '_2' suffixes
    - A separation column is added to the output table showing the distance between matched pairs
    """
    start_time = time.time()
    

    # Check for empty tables
    if len(t1) == 0:
        raise ValueError("First table is empty")
    
    if len(t2) == 0:
        raise ValueError("Second table is empty")
    
    # Check for required columns
    for table_name, table, cols in [("First", t1, t1_coord_cols), ("Second", t2, t2_coord_cols)]:
        for col in cols:
            if col not in table.colnames:
                raise ValueError(f"{table_name} table is missing required column: {col}")
    
    
    # Create SkyCoord objects
    try:
        
        # Suppress warnings about NaN values in coordinates
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Extract coordinate arrays, handling potential issues
            t1_ra = np.array(t1[t1_coord_cols[0]], dtype=float)
            t1_dec = np.array(t1[t1_coord_cols[1]], dtype=float)
            t2_ra = np.array(t2[t2_coord_cols[0]], dtype=float)
            t2_dec = np.array(t2[t2_coord_cols[1]], dtype=float)
            
            # Check for NaN values
            t1_nan_mask = np.isnan(t1_ra) | np.isnan(t1_dec)
            t2_nan_mask = np.isnan(t2_ra) | np.isnan(t2_dec)
            
            t1_nan_count = np.sum(t1_nan_mask)
            t2_nan_count = np.sum(t2_nan_mask)
            
            if t1_nan_count > 0:
                print(f"⚠ First table contains {t1_nan_count} rows with NaN coordinates")
            
            if t2_nan_count > 0:
                print(f"⚠ Second table contains {t2_nan_count} rows with NaN coordinates")
            
            # Create SkyCoord objects
            if t1_nan_count > 0:
                # Handle NaN values by creating a masked array
                valid_t1_mask = ~t1_nan_mask
                t1_coords = SkyCoord(
                    ra=t1_ra[valid_t1_mask]*u.deg, 
                    dec=t1_dec[valid_t1_mask]*u.deg, 
                    frame='icrs'
                )
            else:
                t1_coords = SkyCoord(ra=t1_ra*u.deg, dec=t1_dec*u.deg, frame='icrs')
            
            if t2_nan_count > 0:
                # Handle NaN values by creating a masked array
                valid_t2_mask = ~t2_nan_mask
                t2_coords = SkyCoord(
                    ra=t2_ra[valid_t2_mask]*u.deg, 
                    dec=t2_dec[valid_t2_mask]*u.deg, 
                    frame='icrs'
                )
            else:
                t2_coords = SkyCoord(ra=t2_ra*u.deg, dec=t2_dec*u.deg, frame='icrs')
            
            # Use valid data for t1 and t2
            t1_subset = t1[~t1_nan_mask] if t1_nan_count > 0 else t1
            t2_subset = t2[~t2_nan_mask] if t2_nan_count > 0 else t2
    
    except Exception as e:
        raise
    
    # Perform coordinate matching
    match_start = time.time()
    
    try:
        idx, sep2d, _ = t2_coords.match_to_catalog_sky(t1_coords)
        match_time = time.time() - match_start
    except Exception as e:
        raise
    
    # Apply separation filter
    matched_mask = sep2d < sep * u.arcsec
    match_count = np.sum(matched_mask)
    
    
    if match_count == 0:
        
        # Return empty table with appropriate structure
        merged_cols = t2.colnames + t1.colnames + [f"{match_prefix}sep"]
        merged_table = Table(names=merged_cols, dtype=[t2.dtype[col] for col in t2.colnames] + 
                                                       [t1.dtype[col] for col in t1.colnames] + 
                                                       [float])
        return merged_table
    
    # Get matched rows from both tables
    if t2_nan_count > 0:
        df_matched = t2_subset[matched_mask]
    else:
        df_matched = t2[matched_mask]
        
    tab_selection_matched = t1_subset[idx[matched_mask]]
    
    # Combine matched rows
    try:
        # Resolve column conflicts
        t1_cols = set(tab_selection_matched.colnames)
        t2_cols = set(df_matched.colnames)
        common_cols = t1_cols.intersection(t2_cols)
        
        if common_cols:
            
            # Rename columns in second table to avoid conflicts
            for col in common_cols:
                if col not in t1_coord_cols and col not in t2_coord_cols:
                    df_matched.rename_column(col, f"{col}_2")
        
        # Stack the tables horizontally
        merged_table = hstack([df_matched, tab_selection_matched])
        
        # Add match quality information
        merged_table[f"{match_prefix}sep"] = sep2d[matched_mask].to(u.arcsec)
        
        # Add metadata if requested
        if add_metadata:
            merged_table.meta["MERGETYP"] = "INNER"
            merged_table.meta["MERGESEP"] = sep
            merged_table.meta["MERGECNT"] = match_count
            merged_table.meta["MERGEDAT"] = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())
        
        # Final report
        elapsed = time.time() - start_time
        print(f"✓ Inner merge completed in {elapsed:.2f} seconds")
        print(f"  ↳ Result table: {len(merged_table)} rows, {len(merged_table.colnames)} columns")
        
        return merged_table
        
    except Exception as e:
        raise

# Suppress minor warnings from matplotlib or astroquery (if imported)
warnings.filterwarnings("ignore", category=UserWarning)

# Set Gaia query result limit to a high value
Gaia.ROW_LIMIT = 999999

# Define sky region (RA and DEC boundaries in degrees)
ra = [-100, 100]
dec = [-10, 2]

def main():
    print("Reading S-PLUS DR4 data...")
    
    # Connect to the S-PLUS data API (public DR4)
    conn = splusdata.Core()
    
    print("Querying S-PLUS DR4 data...")
    
    # Query S-PLUS DR4 for sources in the defined RA/DEC range
    # Filters: r-band magnitude < 21, error < 0.03, and CLASS_STAR > 0.9 (likely stars)
    splus_tab = conn.query(
        f"""
        SELECT 
            det.id, det.ra, det.dec, det.field, det.B, det.A, det.class_star, det.theta, det.fwhm,

            -- Photometric errors for all bands
            r.e_r_pstotal, g.e_g_pstotal, i.e_i_pstotal, z.e_z_pstotal, u.e_u_pstotal,
            j0378.e_j0378_pstotal, j0395.e_j0395_pstotal, j0410.e_j0410_pstotal,
            j0430.e_j0430_pstotal, j0515.e_j0515_pstotal, j0660.e_j0660_pstotal, j0861.e_j0861_pstotal,

            -- Fluxes/magnitudes for all bands
            r.r_pstotal, g.g_pstotal, i.i_pstotal, z.z_pstotal, u.u_pstotal,
            j0378.j0378_pstotal, j0395.j0395_pstotal, j0410.j0410_pstotal, 
            j0430.j0430_pstotal, j0515.j0515_pstotal, j0660.j0660_pstotal, j0861.j0861_pstotal

        FROM "dr4_dual"."dr4_dual_detection" AS det 
        
        -- Join with photometric tables for each filter
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_r" AS r ON r.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_g" AS g ON g.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_i" AS i ON i.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_z" AS z ON z.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_u" AS u ON u.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0378" AS j0378 ON j0378.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0395" AS j0395 ON j0395.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0410" AS j0410 ON j0410.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0430" AS j0430 ON j0430.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0515" AS j0515 ON j0515.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0660" AS j0660 ON j0660.id = det.id
        LEFT OUTER JOIN "dr4_dual"."dr4_dual_j0861" AS j0861 ON j0861.id = det.id
        
        WHERE det.ra BETWEEN {ra[0]} AND {ra[1]}
        AND det.dec BETWEEN {dec[0]} AND {dec[1]}

        -- Filters for quality and point-like objects
        AND r.r_pstotal < 21
        AND r.e_r_pstotal < 0.03
        AND CLASS_STAR > 0.9
        """,
        publicdata=True
    )


    print("Querying Gaia DR3 data...")

    # Query Gaia DR3 for sources in the same region
    job = Gaia.launch_job_async(
        "SELECT "
        "source.source_id, source.ra, source.dec, "
        "source.parallax, source.parallax_error, "
        "source.phot_g_mean_mag, source.phot_bp_mean_mag, source.phot_rp_mean_mag, "
        "source.pseudocolour, source.phot_g_mean_flux, source.phot_bp_mean_flux, source.phot_rp_mean_flux, "
        "source.teff_gspphot, source.logg_gspphot, "
        "source.phot_g_mean_flux_error, source.phot_bp_mean_flux_error, source.phot_rp_mean_flux_error, "
        "source.bp_rp, source.bp_g, source.g_rp, "
        "source.phot_proc_mode "
        "FROM gaiadr3.gaia_source AS source "
        f"WHERE source.ra BETWEEN {ra[0]} AND {ra[1]} "
        f"AND source.dec BETWEEN {dec[0]} AND {dec[1]}"
    )

    # Get Gaia results and save to CSV
    r = job.get_results()

    print("Cross-matching S-PLUS and Gaia catalogs...")
    # Cross-match S-PLUS and Gaia sources using inner merge

    # Perform inner merge based on coordinates
    splus_gaia = inner_merge_tables(
        splus_tab,
        r,
        t1_coord_cols=["RA", "DEC"],
        t2_coord_cols=["ra", "dec"],
        sep=1,
        add_metadata=True,
    )

    splus_gaia.write("../data/splus_gaia.csv", overwrite=True, format="csv")
    
    return splus_gaia

# Run the main function
if __name__ == "__main__":
    main()