import csv
import re
import os, sys
import math
import numpy as np
import re
import matplotlib.pylab as plt
import xlsxwriter
import tarfile
import requests
import gzip
import io
from astropy.table import Table
#%matplotlib inline
from IPython.display import Image

from collections import OrderedDict
from scipy.signal import argrelextrema
from string import ascii_lowercase
from os.path import expanduser
from os.path import exists
from astropy.io import ascii, fits
#home = expanduser("~") 
from scipy import signal
#Astropy libraries
from astropy.stats import median_absolute_deviation
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord 
from astropy.coordinates import ICRS, Galactic, Angle, Latitude, Longitude  
from astropy.io import ascii
from astropy.modeling import models
from astropy import constants as const
from astropy.io import fits
#Ignore warnings
import warnings 
warnings.filterwarnings('ignore')
import glob
import time 
home = os.path.expanduser('~')  
