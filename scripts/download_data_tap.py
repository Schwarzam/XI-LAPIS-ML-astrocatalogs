import splusdata
from astroquery.gaia import Gaia

Gaia.ROW_LIMIT = 999999

ra = [-50, 50]
dec = [-2, 2]


def main():
    print("reading splus dr4")
    conn = splusdata.Core()
    
    print("fetching for splus dr4")
    splus_tab = conn.query(
        f"""
        SELECT 
        
        det.id, det.ra, det.dec, det.field, det.B, det.A, det.class_star, det.theta, det.fwhm,
        
        r.e_r_pstotal, g.e_g_pstotal, i.e_i_pstotal, z.e_z_pstotal, u.e_u_pstotal,
        j0378.e_j0378_pstotal, j0395.e_j0395_pstotal, j0410.e_j0410_pstotal,
        j0430.e_j0430_pstotal, j0515.e_j0515_pstotal, j0660.e_j0660_pstotal, j0861.e_j0861_pstotal,
        
        r.r_pstotal, g.g_pstotal, i.i_pstotal, z.z_pstotal, u.u_pstotal,
        j0378.j0378_pstotal, j0395.j0395_pstotal, j0410.j0410_pstotal, 
        j0430.j0430_pstotal, j0515.j0515_pstotal, j0660.j0660_pstotal, j0861.j0861_pstotal
        
        
        
        FROM "dr4_dual"."dr4_dual_detection" AS det 
        
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

        AND r.r_pstotal < 21
        AND r.e_r_pstotal < 0.03
        AND CLASS_STAR > 0.9
        
        """,
        publicdata=True
    )
    splus_tab.write()
    
    print("fetching gaia dr3")
    job = Gaia.launch_job_async("select "
        " "
        "source.source_id, source.ra, source.dec, "
        "source.parallax, source.parallax_error, "
        "source.phot_g_mean_mag, source.phot_bp_mean_mag, source.phot_rp_mean_mag, "
        "source.pseudocolour, source.phot_g_mean_flux, source.phot_bp_mean_flux, source.phot_rp_mean_flux, "
        "source.teff_gspphot, "
        "source.phot_g_mean_flux_error, source.phot_bp_mean_flux_error, source.phot_rp_mean_flux_error, "
        "source.bp_rp, source.bp_g, source.g_rp, "
        "source.phot_proc_mode "
        "from gaiadr3.gaia_source as source "
        f"where source.ra between {ra[0]} and {ra[1]} "
        f"and source.dec between {dec[0]} and {dec[1]} "
    )
    
    r = job.get_results()

    return r

if __name__ == "__main__":
    main()