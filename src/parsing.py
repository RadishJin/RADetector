import gzip
import biotite

with gzip.open("raw_dataset/1CLL.cif.gz", "rt", encoding = "utf-8") as f:
    raw = f.read()

print(raw)