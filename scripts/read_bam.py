import pysam as ps
import pandas as pd

bam = ps.AlignmentFile(
    "/home/isaac/dev/sfu/cmpt415/cmpt415_spatial_domain_variability/data/151507/151507_mRNA.bam",
    'rb'
)

print(bam.fetch().__next__())