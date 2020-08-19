import ReadOIFITS as oifits


dir = '/Users/jacques/Work/Targets/IRAS08544-4431/MATISSE/2020/Iter1_OIFITS_CALIBRATED/MERGED/'
files = '*N_*fits'


data = oifits.read(dir, files)
print(data.vis2)
