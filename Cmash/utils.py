import requests
import numpy as np
import ligo.skymap
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.postprocess import find_greedy_credible_levels