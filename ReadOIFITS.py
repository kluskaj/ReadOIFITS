## OIFITS READING MODULE


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def inform(msg):
    print(bcolors.OKBLUE + msg + bcolors.ENDC)


def inform2(msg):
    print(bcolors.OKGREEN + msg + bcolors.ENDC)


def warn(msg):
    print(bcolors.WARNING + msg + bcolors.ENDC)


def log(msg, dir):
    f = open(dir+"log.txt", "a")
    f.write(msg+"\n")
    f.close()


class data:
    def __init__(self, dir='./', files='*fits', inst='all'):
        self.files = files
        self.dir = dir
        self.inst = inst

    def read(self):
        inform2('Reading from {}{}'.format(self.dir, self.files) )
        inform2('Instrument is {}'.format(self.inst) )
        if inst=='MATISSE':
            readfilesMATISSE(self)

    def readfilesMATISSE(self):
        dir = self.dir
        files = self.files
        listOfFiles = os.listdir(dir)
        i = 0
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, files):
                i += 1
                inform ('Reading '+entry+'...')
                if i == 1:
                    data = readMATISSE(dir+entry)
                else:
                    datatmp = readMATISSE(dir+entry)
                    # Appending all the stuff together
                    # Starting with u coordinates
                    ut, u1t, u2t, u3t = datatmp['u']
                    u, u1, u2, u3 = data['u']
                    u = np.append(u, ut)
                    u1 = np.append(u1, u1t)
                    u2 = np.append(u2, u2t)
                    u3 = np.append(u3, u3t)
                    data['u'] = (u, u1, u2, u3)
                    # v coordinates
                    vt, v1t, v2t, v3t = datatmp['v']
                    v, v1, v2, v3 = data['v']
                    v = np.append(v, vt)
                    v1 = np.append(v1, v1t)
                    v2 = np.append(v2, v2t)
                    v3 = np.append(v3, v3t)
                    data['v'] = (v, v1, v2, v3)
                    # wavelength tables
                    wavet, wavecpt = datatmp['wave']
                    wave, wavecp = data['wave']
                    wave = np.append(wave, wavet)
                    wavecp = np.append(wavecp, wavecpt)
                    data['wave'] = (wave, wavecp)
                    # Visibility squared
                    v2t, v2et = datatmp['v2']
                    v2, v2e = data['v2']
                    v2 = np.append(v2, v2t)
                    v2e = np.append(v2e, v2et)
                    data['v2'] = (v2, v2e)
                    # closure phases
                    cpt, cpet = datatmp['cp']
                    cp, cpe = data['cp']
                    cp = np.append(cp, cpt)
                    cpe = np.append(cpe, cpet)
                    data['cp'] = (cp, cpe)


class OITARGET:
    def __init__(self, target_id=[], target=[]):
        self.target_id = target_id
        self.target = target

    def addtarget(self, target, target_id):
        self.target.extend(target)
        self.target_id.extend(target_id)

    def printtarget(self):
        for t, i in zip(self.target, self.target_id):
            inform('Target #{} is {}'.format(i, t))

    def givetarget(self):
        return self.target, self.target_id

class OIARRAY:
    def __init__(self, arrname='UNKNOWN', tel_name=[], sta_name=[], sta_index=[], diameter=[0]):
        self.arrname = arrname
        self.tel_name = tel_name
        self.sta_name = sta_name
        self.sta_index = sta_index
        self.diameter = diameter

    def addarray(self, arrname, tel_name, sta_name, sta_index, diameter):
        self.arrname.extend(arrname)
        self.tel_name.extend(tel_name)
        self.sta_name.extend(sta_name)
        self.sta_index.extend(sta_index)
        self.diameter.extend(diameter)


class OIWAVE:
    def __init__(self, insname, effwave=[], effband=[]):
        self.insname = insname
        self.effwave = effwave
        self.effband = effband

    def addwave (self, insname, effwave, effband):
        self.effwave.extend(effwave)
        self.effband.extend(effband)
        self.insname.extend(insname)

class OIVIS2:
    def __init__(self, arrname, insname, dateobs=0, mjd=[], vis2data=[], vis2err=[], ucoord=[], vcoord=[], flag=[]):
        self.arrname = arrname
        self.insname = insname
        self.dateobs = dateobs
        self.mjd = mjd
        self.vis2data = vis2data
        self.vis2err = vis2err
        self.ucoord = ucoord
        self.vcoord = vcoord
        self.flag = flag

class OIVIS:
    def __init__(self, arrname, insname, amptype='UNKNOWN', phitype='UNKNOWN', dateobs=0, mjd=[], visamp=[], visamperr=[], visphi=[], visphierr=[], ucoord=[], vcoord=[], flag=[]):
        self.arrname = arrname
        self.insname = insname
        self.dateobs = dateobs
        self.amptype = amptype
        self.phitype = phitype
        self.mjd = mjd
        self.visamp = visamp
        self.visphi = visphi
        self.visamperr = visamperr
        self.visphierr = visphierr
        self.flag = flag

class OIT3:
    def __init__(self, arrname, insname, dateobs=0, mjd=[], t3amp=[], t3amperr=[], t3phi=[], t3phierr=[], u1coord=[], v1coord=[], u2coord=[], v2coord=[], flag=[]):
        self.arrnname = arrname
        self.insname = insname
        self.dateobs = dateobs
        self.mjd = mjd
        self.t3amp = t3amp
        self.t3amperr = t3amperr
        self.t3phi = t3phi
        self.t3phierr = t3phierr
        self.u1coord = u1coord
        self.v1coord = v1coord
        self.u2coord = u2coord
        self.v2coord = v2coord
        self.flag = flag

class OIFLUX:
    def __init__(self, insname, arrname, dateobs=0, mjd=[], fluxdata=[], fluxerr=[], flag=[]):
        self.insname = insname
        self.arrname = arrname
        self.dateobs = dateobs
        self.mjf = mjd
        self.fluxdata = fluxdata
        self.fluxerr = fluxerr
        self.flag = flag
