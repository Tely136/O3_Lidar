import earthaccess
# import os

earthaccess.login(strategy="netrc")

results = earthaccess.search_data(
    short_name='TEMPO_O3PROF_L2',
    version='V04',
    bounding_box=(-74.4, 40.5, -73.5, 41),
    temporal=("2024-10-30", "2024-10-31"),
    count=-1
)

print(results)

files = earthaccess.download(results, "C:/Users/tely1/PhD/o3_processed_data/tempo/")
print(f"Downloaded files: {files}")
