import numpy as np
import warnings
import glob
from obspy import Stream, Trace, UTCDateTime

#
# last update: 25 Jun 2025
#


def __s4(value):
    return -(value & 0b1000) | (value & 0b0111)

def __read_win32_1(filename, chdict):
    output = {}
    srates = {}

    with open(filename, "rb") as fpin:
        fpin.seek(0, 2) # go to end
        sz = fpin.tell() # current position (file length)
        fpin.seek(0) # go to top
        leng = 0
        status0 = 0
        start = 0
        fpin.read(4) # for win32
        leng += 4 # for win32
        while leng < sz:
            # truelen = int.from_bytes(fpin.read(4), 'big')
            # leng += 4
            # if truelen == 0:
            #     break
            fpin.read(1)
            yyyy = "20%02x" % ord(fpin.read(1))
            mm = "%02x" % ord(fpin.read(1))
            dd = "%02x" % ord(fpin.read(1))
            hh = "%02x" % ord(fpin.read(1))
            mi = "%02x" % ord(fpin.read(1))
            ss = "%02x" % ord(fpin.read(1))
            xx = "%02x" % ord(fpin.read(1))
            leng += 8
            date = UTCDateTime(int(yyyy), int(mm), int(dd), int(hh), int(mi), int(ss), int(int(xx)*1e4))

            fpin.read(4) # frame length
            blocklength = int.from_bytes(fpin.read(4), 'big')
            leng += 8

            chlist = [] # to check the duplication
            if start == 0:
                start = date
            if status0 == 0:
                sdata = None
            cloc = 0
            #while leng < truelen:
            while cloc < blocklength:
                fpin.read(1) #oid = fpin.read(1)
                fpin.read(1) #nid = fpin.read(1)
                chnum = ((fpin.read(2)).hex()).lower()
                buff = fpin.read(2)
                datawide = int('%x' % (ord(buff[0:1]) >> 4))
                nsample = int.from_bytes(buff, byteorder='big') & 0b111111111111 #12bit    
                cloc += 6
                data1 = int.from_bytes(fpin.read(4), 'big', signed=True)
                cloc += 4

                # check the duplication
                flag_duplication = chnum in chlist

                if(not flag_duplication):
                    chlist.append(chnum)
                    if chnum in output:
                        output[chnum].append(data1)
                    else:
                        output[chnum] = [data1, ]
                        srates[chnum] = nsample
                len1m = datawide * (nsample-1)
                if datawide == 0:
                    len1m = nsample//2
                sdata = fpin.read(len1m)
                cloc += len1m
                leng += 6+4+len1m

                if len(sdata) < len1m:
                    fpin.seek(-(len1m - len(sdata)), 1)
                    sdata += fpin.read(len1m - len(sdata))
                    msg = "This shouldn't happen, it's weird..."
                    warnings.warn(msg)

                if(flag_duplication):
                    continue

                if datawide == 0:
                    idata2 = data1
                    for i in range(len1m):
                        idata2 = output[chnum][-1] + (int.from_bytes(sdata[i:i + 1], 'big', signed=True) >> 4)
                        output[chnum].append(idata2)
                        if i == (len1m -1):
                            break
                        idata2 += __s4(int.from_bytes(sdata[i:i+1],'big') & 0b1111)
                        output[chnum].append(idata2)
                elif datawide == 1:
                    for i in range((len1m // datawide)):
                        idata2 = output[chnum][-1] + int.from_bytes(sdata[i:i + 1], 'big', signed=True)
                        output[chnum].append(idata2)
                elif datawide == 2:
                    for i in range((len1m // datawide)):
                        idata2 = output[chnum][-1] +\
                                int.from_bytes(sdata[i*2:(i+1)*2], 'big', signed=True)
                        output[chnum].append(idata2)
                elif datawide == 3:
                    for i in range((len1m // datawide)):
                        idata2 = output[chnum][-1] +\
                            int.from_bytes(sdata[i*3:(i+1)*3], 'big', signed=True)
                        output[chnum].append(idata2)
                elif datawide == 4:
                    for i in range((len1m // datawide)):
                        idata2 = output[chnum][-1] +\
                                int.from_bytes(sdata[i*4:(i+1)*4], 'big', signed=True)
                        output[chnum].append(idata2)
                else:
                    msg = "DATAWIDE is %s " % datawide + \
                        "but only values of 0, 1, 2, 3 or 4 are supported."
                    raise NotImplementedError(msg)
    traces = []
    for i in output.keys():
        if(chdict is None):
            t = Trace(data=np.array(output[i]))
        else:
            if(i not in chdict.keys()):
                continue
            t = Trace(data=np.array(output[i])*chdict[i][1])
        t.stats.channel = str(i).lower()
        t.stats.station = chdict[str(i).lower()][0]
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
            #print(line)
            el = decoded_line.split()
            chnum = el[0]
            adstep = float(el[12])
            sens = float(el[7])
            amp = float(el[11])
            if(sens != 0):
                fac = adstep/(sens*pow(10.0, amp/20.0))
            else:
                fac = 0.0
            if(float(el[9]) == 0):
                f0 = 0
            else:
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

def read_win32(*winfiles, chtbl=None):
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
        st += __read_win32_1(filename, chdict)
    #print('Complete!')
    st.merge(method=1)
    if(chtbl is None):
        return st
    else:
        return st, chdict


