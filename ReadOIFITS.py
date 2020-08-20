import astropy.io.fits as fits
import os
import fnmatch
import numpy as np

# OIFITS READING MODULE


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def header(msg):
    print(bcolors.HEADER + msg + bcolors.ENDC)


def bold(msg):
    print(bcolors.BOLD + msg + bcolors.ENDC)


def underline(msg):
    print(bcolors.UNDERLINE + msg + bcolors.ENDC)


def inform(msg):
    print(bcolors.OKBLUE + msg + bcolors.ENDC)


def inform2(msg):
    print(bcolors.OKGREEN + msg + bcolors.ENDC)


def warn(msg):
    print(bcolors.WARNING + msg + bcolors.ENDC)


def fail(msg):
    print(bcolors.FAIL + msg + bcolors.ENDC)


def log(msg, dir):
    f = open(dir+"log.txt", "a")
    f.write(msg+"\n")
    f.close()


def read(dir, files):
    dataset = data(dir, files)
    return dataset


def flatten(L):
    for l in L:
        if isinstance(l, list):
            yield from flatten(l)
        else:
            yield l


class data:
    def __init__(self, dir='./', files='*fits'):
        self.files = files
        self.dir = dir
        self.target = []  # OITARGET()
        self.wave = []  # OIWAVE()
        self.vis2 = []  # OIVIS2()
        self.t3 = []  # OIT3()
        self.vis = []  # OIVIS()
        self.array = []  # OIARRAY()
        self.flux = []  # OIFLUX()
        self.read()
        self.associateWave()
        self.associateFreq()
        self.extendMJD()
        header('Success! \o/')

    def extendMJD(self):
        inform('Assigning mjd...')
        # OIVIS
        for i in np.arange(len(self.vis)):
            mjd = []
            mjd0 = self.vis[i].mjd
            effwave = self.vis[i].effwave
            for j in np.arange(len(mjd0)):
                mjd.append(np.full(len(effwave[j]), mjd0[j]))
            self.vis[i].mjd = np.array(mjd)
        # OIVIS2
        for i in np.arange(len(self.vis2)):
            mjd = []
            mjd0 = self.vis2[i].mjd
            effwave = self.vis2[i].effwave
            for j in np.arange(len(mjd0)):
                mjd.append(np.full(len(effwave[j]), mjd0[j]))
            self.vis2[i].mjd = np.array(mjd)
        # OIT3
        for i in np.arange(len(self.t3)):
            mjd = []
            mjd0 = self.t3[i].mjd
            effwave = self.t3[i].effwave
            for j in np.arange(len(mjd0)):
                mjd.append(np.full(len(effwave[j]), mjd0[j]))
            self.t3[i].mjd = np.array(mjd)

    def giveV2(self, removeflagged=True):
        if self.vis2 == []:
            fail('There is no V2 data in the files')
        else:
            V2, V2e, u, v, lam, mjd = [], [], [], [], [], []
            for data in self.vis2:
                flag = np.logical_not(data.flag)
                V2i = data.vis2data
                V2erri = data.vis2err
                ui = data.uf
                vi = data.vf
                mjdi = data.mjd
                lami = data.effwave

                V2.append(V2i[flag])
                V2e.append(V2erri[flag])
                u.append(ui[flag])
                v.append(vi[flag])
                mjd.append(mjdi[flag])
                lam.append(lami[flag])
        return V2, V2e, u, v, lam, mjd

    def givedataJK(self):
        dataJK = {}
        # OIVIS
        u, v, wave, visamp, visamperr, visphi, visphierr = [], [], [], [], [], [], []
        for i in np.arange(len(self.vis)):
            # fetching
            ui = self.vis[i].uf
            vi = self.vis[i].vf
            wavei = self.vis[i].effwave
            visampi = self.vis[i].visamp
            visamperri = self.vis[i].visamperr
            visphii = self.vis[i].visphi
            visphierri = self.vis[i].visphierr
            # formatting
            ui = flatten(ui)
            vi = flatten(vi)
            wavei = flatten(wavei)
            visampi = flatten(visampi)
            visamperri = flatten(visamperri)
            visphii = flatten(visphii)
            visphierri = flatten(visphierri)
            # loading
            u.extend(ui)
            v.extend(vi)
            wave.extend(wavei)
            visamp.extend(visampi)
            visamperr.extend(visamperri)
            visphi.extend(visphii)
            visphierr.extend(visphierri)
        # flattening and np.arraying
        u = np.array(list(flatten(u)))
        v = np.array(list(flatten(v)))
        wave = np.array(list(flatten(wave)))
        visamp = np.array(list(flatten(visamp)))
        visphi = np.array(list(flatten(visphi)))
        visamperr = np.array(list(flatten(visamperr)))
        visphierr = np.array(list(flatten(visphierr)))
        # writing in the dictionnary
        dataJK['uvV'] = (u.flatten(), v.flatten())
        dataJK['waveV'] = wave.flatten()
        #vis = {}
        #vis['visamp'], vis['visamperr'], vis['visphi'], vis['visphierr'] = visamp, visamperr, visphi, visphierr
        dataJK['vis'] = (visamp.flatten(), visamperr.flatten(), visphi.flatten(), visphierr.flatten())

        # OIVIS2
        u, v, wave, vis2, vis2err = [], [], [], [], []
        for i in np.arange(len(self.vis2)):
            # fetching
            ui = self.vis2[i].uf
            vi = self.vis2[i].vf
            wavei = self.vis2[i].effwave
            vis2i = self.vis2[i].vis2data
            vis2erri = self.vis2[i].vis2err
            # formatting
            ui = flatten(ui)
            vi = flatten(vi)
            wavei = flatten(wavei)
            vis2i = flatten(vis2i)
            vis2erri = flatten(vis2erri)
            # loading
            u.extend(ui)
            v.extend(vi)
            wave.extend(wavei)
            vis2.extend(vis2i)
            vis2err.extend(vis2erri)
        # flattening and np.arraying
        u = np.array(list(flatten(u)))
        v = np.array(list(flatten(v)))
        wave = np.array(list(flatten(wave)))
        vis2 = np.array(list(flatten(vis2)))
        vis2err = np.array(list(flatten(vis2err)))

        # OIT3
        u1, v1, u2, v2, u3, v3, wavecp, t3phi, t3phierr = [], [], [], [], [], [], [], [], []
        for i in np.arange(len(self.vis2)):
            # fetching
            u1i = self.t3[i].uf1
            v1i = self.t3[i].vf1
            u2i = self.t3[i].uf2
            v2i = self.t3[i].vf2
            u3i, v3i = [], []
            for x1, x2, y1, y2 in zip(u1i, u2i, v1i, v2i):
                u3i.extend(x1+x2)
                v3i.extend(y1+y2)
            wavecpi = self.t3[i].effwave
            t3phii = self.t3[i].t3phi
            t3phierri = self.t3[i].t3phierr
            # formatting
            u1i = flatten(u1i)
            v1i = flatten(v1i)
            u2i = flatten(u2i)
            v2i = flatten(v2i)
            u3i = flatten(u3i)
            v3i = flatten(v3i)
            wavecpi = flatten(wavecpi)
            t3phii = flatten(t3phii)
            t3phierri = flatten(t3phierri)
            # loading
            u1.extend(u1i)
            v1.extend(v1i)
            u2.extend(u2i)
            v2.extend(v2i)
            u3.extend(u3i)
            v3.extend(v3i)
            wavecp.extend(wavecpi)
            t3phi.extend(t3phii)
            t3phierr.extend(t3phierri)
        # flattening and np.arraying
        u1 = np.array(list(flatten(u1)))
        v1 = np.array(list(flatten(v1)))
        u2 = np.array(list(flatten(u2)))
        v2 = np.array(list(flatten(v2)))
        u3 = np.array(list(flatten(u3)))
        v3 = np.array(list(flatten(v3)))
        wavecp = np.array(list(flatten(wavecp)))
        t3phi = np.array(list(flatten(t3phi)))
        t3phierr = np.array(list(flatten(t3phierr)))

        # writing in the dictionnary
        dataJK['u'] = (u.flatten(), u1.flatten(), u2.flatten(), u3.flatten())
        dataJK['v'] = (v.flatten(), v1.flatten(), v2.flatten(), v3.flatten())
        dataJK['wave'] = (wave.flatten(), wavecp.flatten())
        dataJK['v2'] = (vis2.flatten(), vis2err.flatten())
        dataJK['cp'] = (t3phi.flatten(), t3phierr.flatten())

        return dataJK

    def associateFreq(self):
        inform('Assigning spatial frequencies...')
        # OIVIS
        for i in np.arange(len(self.vis)):
            uf, vf = [], []
            u = self.vis[i].ucoord
            v = self.vis[i].vcoord
            effwave = self.vis[i].effwave
            for j in np.arange(len(u)):
                uf.append(u[j]/effwave[j])
                vf.append(v[j]/effwave[j])
            self.vis[i].uf = np.array(uf)
            self.vis[i].vf = np.array(vf)
            self.vis[i].base = np.sqrt(np.array(uf)**2 + np.array(vf)**2)
        # OIVIS2
        for i in np.arange(len(self.vis2)):
            uf, vf = [], []
            u = self.vis2[i].ucoord
            v = self.vis2[i].vcoord
            effwave = self.vis2[i].effwave
            for j in np.arange(len(u)):
                uf.append(u[j]/effwave[j])
                vf.append(v[j]/effwave[j])
            self.vis2[i].uf = np.array(uf)
            self.vis2[i].vf = np.array(vf)
            self.vis2[i].base = np.sqrt(np.array(uf)**2 + np.array(vf)**2)
        # OIT3
        for i in np.arange(len(self.t3)):
            uf1, vf1, uf2, vf2 = [], [], [], []
            u1 = self.t3[i].u1coord
            v1 = self.t3[i].v1coord
            u2 = self.t3[i].u2coord
            v2 = self.t3[i].v2coord
            effwave = self.t3[i].effwave
            for j in np.arange(len(u1)):
                uf1.append(u1[j]/effwave[j])
                vf1.append(v1[j]/effwave[j])
                uf2.append(u2[j]/effwave[j])
                vf2.append(v2[j]/effwave[j])
            self.t3[i].uf1 = uf1
            self.t3[i].vf1 = vf1
            self.t3[i].uf2 = uf2
            self.t3[i].vf2 = vf2

    def associateWave(self):
        inform('Assigning wavelengths...')
        # fetch the wavelengths from OIWAVE
        waveid = {}
        for i in np.arange(len(self.wave)):
            name = self.wave[i].insname
            waveid[name] = i

        # Associate the right waves to OIVIS
        for i in np.arange(len(self.vis)):
            name = self.vis[i].insname
            try:
                id = waveid[name]
                effwave = self.wave[id].effwave
                wave = []
                for j in np.arange(self.vis[i].visamp.shape[0]):
                    wave.append(effwave)
                self.vis[i].effwave = np.array(wave)
            except:
                fail('No wavetable corresponding to {} found...'.format(name))

        # Associate the right waves to OIVIS2
        for i in np.arange(len(self.vis2)):
            name = self.vis2[i].insname
            try:
                id = waveid[name]
                effwave = self.wave[id].effwave
                wave = []
                for j in np.arange(self.vis2[i].vis2data.shape[0]):
                    wave.append(effwave)
                self.vis2[i].effwave = np.array(wave)
            except:
                fail('No wavetable corresponding to {} found...'.format(name))

        # Associate the right waves to OIT3
        for i in np.arange(len(self.t3)):
            name = self.t3[i].insname
            try:
                id = waveid[name]
                effwave = self.wave[id].effwave
                wave = []
                for j in np.arange(self.t3[i].t3phi.shape[0]):
                    wave.append(effwave)
                self.t3[i].effwave = np.array(wave)
            except:
                fail('No wavetable corresponding to {} found...'.format(name))

        # Associate the right waves to OIFLUX
        for i in np.arange(len(self.flux)):
            name = self.flux[i].insname
            try:
                id = waveid[name]
                effwave = self.wave[id].effwave
                wave = []
                for j in np.arange(self.flux[i].fluxdata.shape[0]):
                    wave.append(effwave)
                self.flux[i].effwave = np.array(wave)
            except:
                fail('No wavetable corresponding to {} found...'.format(name))

    def read(self):
        header('Reading from {}{}'.format(self.dir, self.files))
        dir = self.dir
        files = self.files
        listOfFiles = os.listdir(dir)
        i = 0
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, files):
                i += 1
                inform('Reading '+entry+'...')
                self.readfile(entry)

    def readfile(self, file):
        hdul = fits.open(self.dir+file)
        err = False
        i = 0
        while err == False:
            i += 1
            try:
                extname = hdul[i].header['EXTNAME']
                print ('Reading '+extname)
                if extname == 'OI_TARGET':
                    self.readTARGET(hdul[i])
                elif extname == 'OI_ARRAY':
                    self.readARRAY(hdul[i])
                elif extname == 'OI_WAVELENGTH':
                    self.readWAVE(hdul[i])
                elif extname == 'OI_VIS':
                    self.readVIS(hdul[i])
                elif extname == 'OI_VIS2':
                    self.readVIS2(hdul[i])
                elif extname == 'OI_T3':
                    self.readT3(hdul[i])
                elif extname == 'OI_FLUX':
                    self.readFLUX(hdul[i])
            except IndexError:
                err = True

    def readTARGET(self, hd):
        target_id = hd.data['TARGET_ID']
        target = hd.data['TARGET']
        tar = OITARGET(target_id=target_id, target=target)
        self.target.append(tar)

    def readARRAY(self, hd):
        arrname = hd.header['ARRNAME']
        tel = hd.data['TEL_NAME']
        sta = hd.data['STA_NAME']
        staid = hd.data['STA_INDEX']
        diam = hd.data['DIAMETER']
        arr = OIARRAY(arrname=arrname, tel_name=tel, sta_name=sta, sta_index=staid, diameter=diam)
        self.array.append(arr)

    def readVIS2(self, hd):
        insname = hd.header['INSNAME']
        arrname = hd.header['ARRNAME']
        dateobs = hd.header['DATE-OBS']
        targetid = hd.data['TARGET_ID']
        mjd = hd.data['MJD']
        vis2data = hd.data['VIS2DATA']
        vis2err = hd.data['VIS2ERR']
        u = hd.data['UCOORD']
        v = hd.data['VCOORD']
        sta = hd.data['STA_INDEX']
        flag = hd.data['FLAG']
        vis2 = OIVIS2(arrname, insname, dateobs=dateobs, mjd=mjd, vis2data=vis2data, vis2err=vis2err, ucoord=u, vcoord=v, flag=flag)
        self.vis2.append(vis2)

    def readT3(self, hd):
        insname = hd.header['INSNAME']
        arrname = hd.header['ARRNAME']
        dateobs = hd.header['DATE-OBS']
        t3amp = hd.data['T3AMP']
        t3phi = hd.data['T3PHI']
        t3amperr = hd.data['T3AMPERR']
        t3phierr = hd.data['T3PHIERR']
        mjd = hd.data['MJD']
        targetid = hd.data['TARGET_ID']
        u1 = hd.data['U1COORD']
        u2 = hd.data['U2COORD']
        v1 = hd.data['V1COORD']
        v2 = hd.data['V2COORD']
        staid = hd.data['STA_INDEX']
        flag = hd.data['FLAG']
        T3 = OIT3(arrname, insname, dateobs=dateobs, mjd=mjd, t3amp=t3amp, t3amperr=t3amperr, t3phi=t3phi, t3phierr=t3phierr, u1coord=u1, v1coord=v1, u2coord=u2, v2coord=v2, flag=flag)
        self.t3.append(T3)

    def readVIS(self, hd):
        insname = hd.header['INSNAME']
        arrname = hd.header['ARRNAME']
        try:
            amptype = hd.header['AMPTYP']
            phitype = hd.header['PHITYP']
        except KeyError:
            amptype = 'unknown'
            phitype = 'unknown'
        dateobs = hd.header['DATE-OBS']
        mjd = hd.data['MJD']
        targetid = hd.data['TARGET_ID']
        visamp = hd.data['VISAMP']
        visamperr = hd.data['VISAMPERR']
        visphi = hd.data['VISPHI']
        visphierr = hd.data['VISPHIERR']
        u = hd.data['UCOORD']
        v = hd.data['VCOORD']
        staid = hd.data['STA_INDEX']
        flag = hd.data['FLAG']
        VIS = OIVIS(arrname, insname, amptype=amptype, phitype=phitype, dateobs=dateobs, mjd=mjd, visamp=visamp, visamperr=visamperr, visphi=visphi, visphierr=visphierr, ucoord=u, vcoord=v, flag=flag)
        self.vis.append(VIS)

    def readWAVE(self, hd):
        insname = hd.header['INSNAME']
        effwave = np.array(hd.data['EFF_WAVE'])
        effband = np.array(hd.data['EFF_BAND'])
        wave0 = OIWAVE(insname, effwave=effwave, effband=effband)
        self.wave.append(wave0)
        # print(self.wave, self.wave[0].effwave)

    def readFLUX(self, hd):
        dateobs = hd.header['DATE-OBS']
        insname = hd.header['INSNAME']
        arrname = hd.header['ARRNAME']
        calstat = hd.header['CALSTAT']
        targetid = hd.data['TARGET_ID']
        mjd = hd.data['MJD']
        try:
            flux = hd.data['FLUXDATA']
        except KeyError:
            flux = hd.data['FLUX']  #  for GRAVITY
        fluxerr = hd.data['FLUXERR']
        staid = hd.data['STA_INDEX']
        flag = hd.data['FLAG']
        fl = OIFLUX(insname, arrname, calstat=calstat, dateobs=dateobs, mjd=mjd, fluxdata=flux, fluxerr=fluxerr, flag=flag)
        self.flux.append(fl)




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

    def givetargetid(self):
        return self.target, self.target_id

    def givetarget(self):
        return self.target

    def giveid(self):
        return self.target_id

    def givetheid(self, target):
        id = []
        for t, i in zip(self.target, self.target_id):
            if t==target:
                id.extend(i)
        return id

    def givethetarget(self, id):
        tar = []
        for t, i in zip(self.target, self.target_id):
            if i==id:
                tar.extend(t)
        return tar


class OIARRAY:
    def __init__(self, arrname='UNKNOWN', tel_name=[], sta_name=[], sta_index=[], diameter=[]):
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
        self.ucoord = ucoord
        self.vcoord = vcoord
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
    def __init__(self, insname, arrname, calstat='unknown', dateobs=0, mjd=[], fluxdata=[], fluxerr=[], flag=[]):
        self.insname = insname
        self.arrname = arrname
        self.dateobs = dateobs
        self.mjf = mjd
        self.fluxdata = fluxdata
        self.fluxerr = fluxerr
        self.flag = flag
        self.calstat = calstat
