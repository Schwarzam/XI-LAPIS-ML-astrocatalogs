import lsdb
from astroquery.gaia import Gaia
from dask.distributed import Client

# Set Gaia row limit to ensure large query capacity (if used directly from astroquery)
Gaia.ROW_LIMIT = 999999

def main():
    print("Reading S-PLUS DR4 from LSDB...")

    # Load S-PLUS DR4 HATS catalog from LSDB (Dual mode) using `read_hats`
    # Optionally preload margins for quick local spatial indexing
    splus_dr4 = lsdb.read_hats(
        path='https://splus.cloud/HIPS/catalogs/dr4/dual', 
        margin_cache='https://splus.cloud/HIPS/catalogs/dr4/dual_2arcsec',
        columns=[
            "ID", "RA", "DEC", "FWHM",
            "r_PStotal", "e_r_PStotal",
            "g_PStotal", "e_g_PStotal",
            "i_PStotal", "e_i_PStotal",
            "z_PStotal", "e_z_PStotal",
            "u_PStotal", "e_u_PStotal",
            "J0378_PStotal", "e_J0378_PStotal",
            "J0395_PStotal", "e_J0395_PStotal",
            "J0410_PStotal", "e_J0410_PStotal",
            "J0430_PStotal", "e_J0430_PStotal",
            "J0515_PStotal", "e_J0515_PStotal",
            "J0660_PStotal", "e_J0660_PStotal",
            "J0861_PStotal", "e_J0861_PStotal",
            "CLASS_STAR"
        ],
        filters=[
            ("r_PStotal", "<", 21),
            ("e_r_PStotal", "<", 0.03),
            ("CLASS_STAR", ">", 0.9),
        ]
    )

    print("Reading Gaia DR3 from LSDB...")
    
    # Load Gaia DR3 catalog from LSDB via HATS structure
    gaia = lsdb.read_hats(
        path='https://data.lsdb.io/hats/gaia_dr3/gaia', 
        margin_cache='https://data.lsdb.io/hats/gaia_dr3/gaia_10arcs',
        columns=["teff_gsphot"]
    )

    print("Cross-matching S-PLUS and Gaia catalogs...")
    
    # Crossmatch S-PLUS and Gaia sources using internal spatial matching
    splus_gaia = splus_dr4.crossmatch(
        gaia, 
        suffixes=("_dual", "_gaia")  # To avoid column name collisions
    )

    print("Filtering crossmatch by box search...")
    
    # Apply a rectangular region-of-interest selection: RA from -50 to 50, Dec from -2 to 2
    table = splus_gaia.box_search(
        [-50, 50],
        [-2, 2]
    ).compute()  # Convert lazy dask dataframe into in-memory pandas dataframe

    # Save resulting matched table to CSV
    table.to_csv("../data/splus_gaia.csv", index=False)
    print("Done. Result saved to splus_gaia.csv.")

if __name__ == '__main__':
    # Initialize Dask cluster with 10 workers and 8GB RAM limit per worker
    client = Client(n_workers=10, memory_limit="8GB")
    
    # Print dashboard link to monitor Dask progress
    dashboard_link = client.dashboard_link
    print(f"Check the status of the Dask cluster at {dashboard_link}")
    
    main()