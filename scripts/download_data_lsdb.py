import lsdb
from astroquery.gaia import Gaia
from dask.distributed import Client

Gaia.ROW_LIMIT = 999999

def main():
    print("reading splus dr4")
    splus_dr4 = lsdb.read_hats(
        path='https://splus.cloud/HIPS/catalogs/dr4/dual', 
        margin_cache='https://splus.cloud/HIPS/catalogs/dr4/dual_2arcsec'
    )

    splus_dr4 = lsdb.read_hats(
        path='https://splus.cloud/HIPS/catalogs/dr4/dual', 
        margin_cache='https://splus.cloud/HIPS/catalogs/dr4/dual_2arcsec',
        columns=["ID", "RA", "DEC", "FWHM",
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
                "CLASS_STAR", #"EBV_SCH"
            ],
    )

    print("reading gaia dr3")
    gaia = lsdb.read_hats(
        path = 'https://data.lsdb.io/hats/gaia_dr3/gaia', 
        margin_cache='https://data.lsdb.io/hats/gaia_dr3/gaia_10arcs'
    )

    print("cross matching gaia and splus")
    splus_gaia = splus_dr4.crossmatch(
        gaia, 
        suffixes=("_dual", "_gaia"),
    )

    # making a box search 
    table = splus_gaia.box_search(
        [-50, 50],
        [-2, 2]
    ).compute()

    table.to_csv("splus_gaia.csv", index=False)
    print("done")

if '__main__' == __name__:
    # Start Dask client
    client = Client(n_workers=10, memory_limit="8GB")

    dashboard_link = client.dashboard_link
    print(f"Check the status of the Dask cluster at {dashboard_link}")
    
    main()