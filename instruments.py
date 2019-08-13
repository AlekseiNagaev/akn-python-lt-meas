import visa
import re
import time
from math import sqrt,atan2
rm = visa.ResourceManager()
rm.list_resources()

class instr:

    def __init__(obj, str):
        obj.addr=str
        p1 = re.compile('\d{1,3}[.]\d{1,3}[.]\d+[.]\d+')
        p2 = re.compile('COM[0-9]')
        if re.search(p1, obj.addr):
            obj.dev = rm.open_resource('TCPIP::%s::INSTR' % obj.addr)
        elif re.search(p2, obj.addr):
            obj.dev = rm.open_resource('%s' % obj.addr)
        else:
            obj.dev = rm.open_resource('GPIB::%s' % obj.addr)
        if  type(obj) is rbrg:
            obj.model = 'AVS'
        else:
            obj.model = obj.mod()

    def __str__(obj):
        return "Device %s with address %s" % (obj.model, obj.addr)

    def close(obj):
        obj.dev.control_ren(6)
        obj.dev.close()

    def send(obj,s):
        obj.dev.write(s)

    def get(obj):
        x = obj.dev.read()
        return x

    def rst(obj):
        obj.send('*RST')

    def err(obj):
        obj.send('SYST:ERR?')
        s = obj.get()
        print(obj.model + ' ' + s)

    def mod(obj):
        obj.send('*IDN?')
        s = obj.get()
        s = s.split(',')[1]
        #print(s)
        return s

class wfsg(instr):

    r_ch =[0.3169,0.6319,1.259,2.519,5.019]
    devices = ['33120A','33509B','33522A']

    def __init__(obj, str, res = 1e6):
        instr.__init__(obj,str)
        obj.res = res
        obj.bum_per  = 1
        obj.bum_ncyc = 1
        obj.bum_phas = 0
        if obj.model is wfsg.devices[0]:
            obj.gg = 0
        elif obj.model in wfsg.devices[1:2]:
            obj.gg = 1
        else:
            raise ValueError('Instrument %s is probably not a supported wfsg' % obj.model)

    def trg(obj,tc):
        if tc:
            obj.send('*TRG')

    def set_shp(obj,sh):
        if obj.gg:
            obj.send('FUNC:SHAP %s' % sh)
        else:
            obj.send('FUNC %s' % sh)

    def apply(obj, shp, frq, amp, ofs):
        obj.send('APPL:%s %f,%f,%f' % (shp,frq,amp,ofs))

    def set_dcl(obj, dcl = 50):
        if obj.gg:
            obj.send('PULS:DCYC %d' % dcl)
        else:
            obj.send('FUNC:SQU:DCYC %d' % dcl)

    def set_frq(obj, frq = 1e3):
        obj.send('FREQ %f' % frq)

    def set_amp(obj, amp = 1):
        if type(amp) is str:
            obj.send('VOLT %s' % amp)
        else:
            obj.send('VOLT %f' % amp)

    def set_ofs(obj, ofs = 0):
        obj.send('VOLT:OFFS %f' % ofs)

    def get_dc(obj):
        obj.set_shp('DC')
        obj.send('VOLT:OFFS?')
        return float(obj.get())

    def set_dc(obj, vol, p=0.1, ivb='vol'):
        if ivb is 'vol':
            obj.set_out_load('DEF')
        elif ivb is 'cur':
            obj.set_out_load()
        vol1 = obj.get_dc()
        obj.set_out(1)
        obj.set_ofs(vol)
        if obj.model is wfsg.devices[0]:
            if any(abs(vol1) <= r <= abs(vol) for r in wfsg.r_ch) or any(abs(vol1) >= r >= abs(vol) for r in wfsg.r_ch):
                time.sleep(1)
        elif obj.model in wfsg.devices[1:2]:
            if (abs(vol1) <= wfsg.r_ch[0] <= abs(vol)) or (abs(vol1) >= wfsg.r_ch[0] >= abs(vol)):
                time.sleep(1)
        time.sleep(p)

    def set_gate(obj,vol,stp,p=0.1,ivb='vol'):
        vol1 = obj.get_dc()
        if vol is not vol1:
            n = abs(round((vol-vol1)/stp))
            for i in range(n):
                vol2 = vol1+(vol-vol1)*(i+1)/n
                obj.set_dc(vol2,0.05,ivb)
            obj.set_dc(vol,0,ivb)
            print('Vg = %f' % vol)
            time.sleep(p)

    def set_out(obj,fl):
        if fl:
            obj.send('OUTP 1')
        else:
            obj.send('OUTP 0')

    def set_out_load(obj,l=1):
        if obj.model == '33509B':
            if l:
                obj.send('OUTP:LOAD INF')
            else:
                obj.send('OUTP:LOAD %s' % l)

    def set_bum_par(obj):
        if not obj.gg:
            obj.send('BM:INT:RATE %f' % 1/obj.bum_per)
            obj.send('BM:NCYC %d' % obj.bum_ncyc)
            obj.send('BM:PHAS %f' % obj.bum_phas)
        else:
            obj.send('BURS:MODE TRIG')
            obj.send('BURS:INT:PER %f' % obj.bum_per)
            obj.send('BURS:NCYC %d' % obj.bum_ncyc)
            obj.send('BURS:PHAS %f' % obj.bum_phas)

    def set_bum(obj,fl):
        if not obj.gg:
            if fl:
                obj.send('BM:STAT ON')
            else:
                obj.send('BM:STAT OFF')
        else:
            if fl:
                obj.send('BURS:STAT ON')
            else:
                obj.send('BURS:STAT OFF')

    def set_pul_par(obj,frq=1e3,dcl=50):
        obj.set_shp('SQU')
        obj.set_frq(frq)
        obj.set_dcl(dcl)
        obj.set_amp('MIN')
        obj.set_ofs(0);
        obj.set_out_load();
        obj.set_bum_par();
        obj.send('TRIG:SOUR BUS')
        obj.set_bum(1);
        obj.set_out(1);

    def set_pul_amp(obj,amp):
        obj.set_amp(amp)
        obj.set_ofs(amp/2)

class mult(instr):

    res = [1e-5, 1e-7, 3e-8]*10
    npl = [0.2, 1, 10]
    rng = [1e-1, 1, 10]
    dev = ['34401A','MODEL2450']

    def __init__(obj,str,num=1):
        instr.__init__(obj,str)
        obj.num = num
        if not (obj.model in mult.dev):
            raise ValueError('Instrument %s is probably not a supported mult' % obj.model)

    def set_zer(obj,fl):
        if fl:
            obj.send('CALC:FUNC NULL')
            obj.send('CALC:STAT ON')
        else:
            obj.send('CALC:STAT OFF')

    def read_vm(obj):
        obj.send('INIT')
        obj.send('FETCH?')
        v = float(obj.get())
        return v

    def set_vm(obj,num=4):
        obj.set_dc_rng(1)
        time.sleep(0.5)
        time.sleep(0.5)
        if num in range(4,7):
            obj.send('VOLT:DC:RES %.5f' % mult.res[num-4])
            obj.send('VOLT:DC:NPLC %f' % mult.npl[num-4])
        else:
            obj.send('VOLT:DC:RES %.5f' % mult.res[0])
            obj.send('VOLT:DC:NPLC %f' % mult.npl[0])
        obj.send('ZERO:AUTO OFF')
        obj.send('TRIG:SOUR IMM')

    def set_dc_rng(obj,rng):
        obj.send('VOLT:DC:RANGe %f' % rng)

    def read_vm_avg(obj, q, fl):
        v0 = 0
        for i in range(obj.num):
            v = obj.read_vm()
            v0 = v0 + v
        v0 = v0/obj.num
        for i in range(len(mult.rng)):
            if (((abs(v0) > mult.rng[i]) and q) or ((abs(v0) < mult.rng[i]) and not q)) and (fl is q):
                if q:
                    obj.set_dc_rng(mult.rng[i+1])
                else:
                    obj.set_dc_rng(mult.rng[i])
                fl = not q
        time.sleep(1)
        yield v0
        yield fl

class rbrg(instr):

    cha = [0, 1, 2, 3, 4, 5, 6, 7];
    exc = [5, 5, 5, 5, 0, 0, 0, 0];
    # Conversion coefficients
    A = 0.9471
    B = 6.653
    def __init__(obj,str,mod=1):
        instr.__init__(obj,str)
        obj.mod = mod

    def set_rmt(obj,fl):
        if fl:
            obj.send('REM 1')
        else:
            obj.send('REM 0')

    def get_brg(obj):
        obj.send('CH?;RAN?;EXC?')
        s = obj.get()
        return s

    def set_brg_0(obj):
        obj.send('CH?');
        s = obj.get()
        obj.send('REFID0;CH0')

    def set_brg(obj,chn=4,rng=4,exc=2):
        if not((chn in range(7)) and (rng in range(7)) and (exc in range(7))):
            raise ValueError('Wrong parameters for the rbrg. Configuration not changed.')

    def brg_read(obj,avg=1): # Finish r2t()
        if obj.mod:
            avg = max(1,round(avg))
            try:
                obj.send('RES %d;RES?' % avg)
                time.sleep(.2*avg)
                R = float(obj.get())
                #T = r2t(R)
            except:
                print('Warning: reading resistance bridge at serial port ' + obj.addr + ' failed.')
                R = 0
                #T = 0
        else:
            R = 0
            #T = 0
        yield R
        #yield T

class dcpw(instr):

    dev = ['E3647A']

    def __init__(obj,str):
        instr.__init__(obj,str)
        if not (obj.model in dcpw.dev):
            raise ValueError('Instrument %s is probably not a supported dcpw' % obj.model)

    def set_amp(obj,vol,p=1):
        obj.send('VOLT %f' % vol)
        obj.send('OUTP ON')
        time.sleep(p)

    def get_amp(obj):
        obj.send('VOLT?')
        return float(obj.get())

    def set_out(obj,num=1):
        obj.send('INST OUTP %d' % num)

    def set_gate(obj,vol,stp,fl=1,p=1):
        obj.set_out()
        vol1 = obj.get_amp()
        if vol is not vol1:
            n = abs(round((vol-vol1)/stp))
            for i in range(n):
                vol2 = vol1+(vol-vol1)*(i+1)/n
                obj.set_amp(vol2,0.05)
            obj.set_amp(vol,0)
            print('Vg = %f' % vol)
            time.sleep(p)

class liao(instr):

    def __init__(obj,str):
        instr.__init__(obj,str)

    def read_x_y(obj):
        obj.send('X.')
        Vx = float(obj.get())
        obj.send('Y.')
        Vy = float(obj.get())
        R = sqrt(Vx**2 + Vy**2)
        ang = atan2(Vy,Vx)
        yield R
        yield ang
