#!/home/nannu1375/miniconda3/envs/knc_signature/bin/python

# =============================================================================
# 06_spatial_analysis.py
#
# Spatial transcriptomics analysis of Visium HD LUAD dataset
# =============================================================================

# Importing packages
import os
from pathlib import Path

import pandas as pd
import numpy as np

import scanpy as sc
import anndata as ad

import matplotlib.pyplot as plt

# Paths to files

from pathlib import Path

project = Path("/media/nannu1375/Backpack/Shankara/KNC")

spatial_dir = project / "public_data" / "Spatial_data" / "Output-suppli_files"

bin8 = spatial_dir / "binned_outputs" / "square_008um"

print(bin8)