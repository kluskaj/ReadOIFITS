import ReadOIFITS as oifits
import inspect
import matplotlib.pyplot as plt
import numpy as np

# MATISSE N#
dir = '/Users/jacques/Work/Targets/IRAS08544-4431/MATISSE/2020/Iter1_OIFITS_CALIBRATED/MERGED/'
files = '*N_*fits'


# MATISSE L#
files = '*LM_*fits'
# data = oifits.read(dir, files)


# PIONIER#
dir = '/Users/jacques/Work/Targets/IRAS08544-4431/data/'
files = '*fits'


# GRAVITY#
dir = '/Users/jacques/Work/Targets/IRAS08544-4431/GRAVITY/DATA/'
files = 'SCI_IRAS08544-4431_A_SINGLE_SCI_VIS_CALIBRATED_1.fits'


# MIDI#
dir = '/Users/jacques/Work/Targets/IRAS08544-4431/MIDI/'
files = 'IRAS08544-4431_1.2015-01-21b.calvisibility.fits'
# data = oifits.read(dir, files)


# SPHERE/SAM#
dir = '/Users/jacques/Work/INSPIRING/data/SPHERE/'
files = 'IW_Car_0016mrg_allcal.oifits'
# data = oifits.read(dir, files)


# NIRC2/SAM#
dir = '/Users/jacques/Work/Targets/ZCMa/SAM/'
files = 'ZCMa_0363mrg.oifits'
# data = oifits.read(dir, files)


# MIRCX#
dir = '/Users/jacques/Work/Targets/CLLac/'
files = '*fits'
# data = oifits.read(dir, files)
data = oifits.read(dir, files)
data.plotV2CP(lines=True)
data.plotV2CP(lines=False)
