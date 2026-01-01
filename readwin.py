import numpy as np
import warnings
import glob
from obspy import Stream, Trace, UTCDateTime

#
# last update: 11 Nov 2024
# EUC-JP encoding is supported
# Read f0 and h
#


def __s4(value):
    return -(value & 0b1000) | (value & 0b0111)

def __read_win_1(filename, chdict):
    output = {}
    srates = {}

    with open(filename, "rb") as fpin:
        fpin.seek(0, 2)
        sz = fpin.tell()
        fpin.seek(0)
        leng = 0
        status0 = 0
        start = 0
        while leng < sz:
            truelen = int.from_bytes(fpin.read(4), 'big')
            leng = 4
            if truelen == 0:
                break
            buff = fpin.read(6)
            leng += 6

            yy = "20%02x" % ord(buff[0:1])
            mm = "%02x" % ord(buff[1:2])
            dd = "%02x" % ord(buff[2:3])
            hh = "%02x" % ord(buff[3:4])
            mi = "%02x" % ord(buff[4:5])
            sec = "%02x" % ord(buff[5:6])

            date = UTCDateTime(int(yy), int(mm), int(dd), int(hh), int(mi), int(sec))

            chlist = [] # to check the duplication
            if start == 0:
                start = date
            if status0 == 0:
                sdata = None
            while leng < truelen:
                buff = fpin.read(4)
                leng += 4
                flag = '%02x' % ord(buff[0:1])
                chanum = '%02x' % ord(buff[1:2])
                chanum = "%02s%02s" % (flag, chanum)
                datawide = int('%x' % (ord(buff[2:3]) >> 4))
                #srate = ord(buff[3:4])
                srate = int(bin(int.from_bytes(buff, byteorder='big'))[-12:], 2)
                xlen = (srate - 1) * datawide
                if datawide == 0:
                    xlen = srate // 2
                    datawide = 0.5

                leng += 4
                idata22 = int.from_bytes(fpin.read(4), 'big', signed=True)

                # check the duplication
                flag_duplication = chanum in chlist

                if(not flag_duplication):
                    chlist.append(chanum)
                    if chanum in output:
                        output[chanum].append(idata22)
                    else:
                        output[chanum] = [idata22, ]
                        srates[chanum] = srate

                sdata = fpin.read(xlen)
                leng += xlen

                if len(sdata) < xlen:
                    fpin.seek(-(xlen - len(sdata)), 1)
                    sdata += fpin.read(xlen - len(sdata))
                    msg = "This shouldn't happen, it's weird..."
                    warnings.warn(msg)

                if(flag_duplication):
                    continue

                if datawide == 0.5:
                    for i in range(xlen):
                        idata2 = output[chanum][-1] + (int.from_bytes(sdata[i:i + 1], 'big', signed=True) >> 4)
                        output[chanum].append(idata2)
                        if i == (xlen -1):
                            break
                        idata2 += __s4(int.from_bytes(sdata[i:i+1],'big') & 0b1111)
                        output[chanum].append(idata2)
                elif datawide == 1:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] + int.from_bytes(sdata[i:i + 1], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 2:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                                int.from_bytes(sdata[i*2:(i+1)*2], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 3:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                            int.from_bytes(sdata[i*3:(i+1)*3], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 4:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                                int.from_bytes(sdata[i*4:(i+1)*4], 'big', signed=True)
                        output[chanum].append(idata2)
                else:
                    msg = "DATAWIDE is %s " % datawide + \
                        "but only values of 0.5, 1, 2, 3 or 4 are supported."
                    raise NotImplementedError(msg)
    traces = []
    for i in output.keys():
        if(chdict is None):
            t = Trace(data=np.array(output[i]))
        else:
            if(i not in chdict.keys()):
                continue
            t = Trace(data=np.array(output[i])*chdict[i][1])
            t.stats.station = chdict[str(i).lower()][0]
        t.stats.channel = str(i).lower()
        t.stats.sampling_rate = float(srates[i])
        t.stats.starttime = start
        traces.append(t)
    return Stream(traces=traces)

def read_chtable(chtbl):
    chdict = {}
    with open(chtbl, 'rb') as f:
        for line in f:
            try:
                decoded_line = line.decode('utf-8')
            except UnicodeDecodeError:
                decoded_line = line.decode('euc-jp')
            if('#' in decoded_line):
                continue
            #el = line.split()
            el = decoded_line.split()
            chnum = el[0]
            adstep = float(el[12])
            sens = float(el[7])
            amp = float(el[11])
            if(sens != 0):
                fac = adstep/(sens*pow(10.0, amp/20.0))
            else:
                fac = 0.0
            f0 = 1/float(el[9])
            damp = float(el[10])
            stcode = el[3]
            comp = el[4]
            stname = stcode + '.' + comp
            if(len(el) >= 15):
                stlon = float(el[14])
                stlat = float(el[13])
                chdict[chnum.lower()] = [stname.lower(), fac, stcode, f0, damp, stlon, stlat]
            else:
                chdict[chnum.lower()] = [stname.lower(), fac, stcode, f0, damp]
    return chdict

def read_win(*winfiles, chtbl=None):
    if(chtbl is None):
        chdict = None
    else:
        chdict = read_chtable(chtbl)
    winfile = []
    st = Stream()
    if(type(winfiles[0]) is list):
        for path in winfiles[0]:
            print(path)
            filenames = sorted(glob.glob(path))
            for filename in filenames:
                winfile.append(filename)
    else:
        for path in winfiles:
            print(path)
            filenames = sorted(glob.glob(path))
            for filename in filenames:
                winfile.append(filename)
    #print('Reading:', end=' ')
    for filename in winfile:
        #print(filename, end=', ')
        st += __read_win_1(filename, chdict)
    #print('Complete!')
    st.merge(method=1)
    if(chtbl is None):
        return st
    else:
        return st, chdict


