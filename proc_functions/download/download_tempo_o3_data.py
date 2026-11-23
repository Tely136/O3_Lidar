import earthaccess
import os

root_dir = ""
data_dir = "/o3_processed_data/tempo/"

earthaccess.login(strategy="netrc")

results = earthaccess.search_data(
    short_name='TEMPO_O3PROF_L2',
    version='V04',
    # bounding_box=(-10, 20, 10, 50),
    temporal=("2026-02", "2026-03"),
    count=1
)

files = earthaccess.download(results, "C:/Users/tely1/PhD/o3_processed_data/tempo/")
print(f"Downloaded files: {files}")
