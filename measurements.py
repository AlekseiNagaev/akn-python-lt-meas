def ivms(dev,opt,ivp):
    """IV measurement program."""

    # Importing modules
    from instruments import wfsg,mult,dcpw,rbrg         # import instrument classes
    import numpy as np                                  # import numpy for array operations

    if 'wgm' in dev:                                    # set main waveform generator
        wgm = wfsg(dev['wgm'])
    else:
        raise ValueError('No main wfsg')

    if 'wgg' in dev:                                    # set gate generator
        wgg = wfsg(dev['wgg'])
    elif 'psg' in dev:
        wgg = dcpw(dev['psg'])
    else:
        raise ValueError('No gate generator')

    if 'mu1' in dev:                                    # set main multimeter
        mu1 = mult(dev['mu1'])
    else:
        raise ValueError('No main mult')

    if 'mu2' in dev:                                    # set second multimeter?
        mu2 = mult(dev['mu2'])

    if 'brg' in dev:                                    # set resistance bridge
        brg = rbrg(dev['brg'][0],dev['brg'][1])

    # General settings
    bre = opt['bre']                                    # Bias resistance
    gnv = opt['gnv']                                    # gain for voltage amplifier
    gnc = opt['gnc']                                    # gain for current amplifier
    ggp = opt['ggp']                                    # gate generator pause
    swl = opt['swl']                                    # switching level for the jump to normal state
    gss = opt['gss']                                    # gate sweep speed
    ivb = opt['ivb']                                    # bias type
    gcs = opt['gcs']                                    # Gate on/off

    # IV settings
    ivs = ivp['ivs']                                    # iv sweep speed/pause
    pol = ivp['pol']                                    # iv polarity
    amp = ivp['amp']                                    # Max amplitude for IV
    stp = ivp['stp']                                    # sweep step for IV
    ggv = ivp['ggv']                                    # gate voltage values np.arange(x,y,st)

    # Setting up arrays
    lgv = len(ggv)                                      # length of gate voltage array
    vol = np.concatenate((np.arange(0,amp,stp),[amp]))  # Voltage array for the wgm
    lev = len(vol)                                      # Length of wgm voltage array
    data = {'vol': np.zeros((lgv,4,lev)),               # Measured voltage array
            'cur': np.zeros((lgv,4,lev)),               # "Measured" current array
            'res': np.zeros(lgv),                       # Calculated resistance array
            'isw': np.zeros(lgv),                       # Switching current array
            'irt': np.zeros(lgv)                        # Retrapping current array
            }
    # Save file
    # Figures
    # Measurement
    try:
        # Read temperature
        # Presetting instruments
        mu1.set_vm(5)
        if ivb is 'cur':
            mu2.set_vm(5)
        wgm.set_gate(0,gss,1,0,ivb)                     # Set wgm to 0V DC carefully

        # Gate cycle
        for i in range(lgv):
            fl1 = 0
            fl2 = 0
            wgg.set_gate(ggv[i],gss,gcs,ggp,ivb)        # Set gate voltage
            mu1.set_zer(1)                              # Set zero of mu1
            if ivb is 'cur':
                mu2.set_zer(1)                          # Set zero of mu2

            # IV cycle
            for j in range(4):
                if j is 1:
                    vol = vol[::-1]                     # Swap direction
                elif j is 2:
                    vol = - vol[::-1]                   # Swap direction and sign
                elif j is 3:
                    vol = vol[::-1]                     # Swap direction
                # Branch cycle
                for k in range(lev):
                    wgm.set_dc(vol[k],ivs,ivb)          # Set voltage on wgm
                    data['vol'][i,j,k],fl1 = mu1.read_vm_avg(j % 2,fl1)
                    # Change sign if polarity is wrong
                    if pol:
                        data['vol'][i,j,k] = - data['vol'][i,j,k]/gnv
                    else:
                        data['vol'][i,j,k] = data['vol'][i,j,k]/gnv

                    # Calculate/measure applied current
                    if ivb is 'cur':
                        data['cur'][i,j,k],fl2 = mu2.read_vm_avg(j % 2,fl2)
                        data['cur'][i,j,k] = data['cur'][i,j,k]/gnc
                    elif ivb is 'vol':
                        data['cur'][i,j,k] = vol[k]/bre

                    # Plot data

                    # Append to file
            mu1.set_zer(0)
            data['res'][i] = data['vol'][i,0,lev-1]/data['cur'][i,0,lev-1] # Calculate resistance
            print('RES = %.3f' % data['res'][i])
            # Final plot
            # Stop measurements

        # Set gate to 0
        wgg.set_gate(0,gss,gcs,0,ivb)
        # Read temperature
        # Save file
    except:
        raise
    finally:
        wgm.set_gate(0,0.05)
        wgm.close()
        wgg.close()
        mu1.close()
        if 'mu2' in dev:
            mu2.close()

    # Switching and retrapping current calculations
    for i in range(lgv):
        try:
            j = np.argwhere(data['vol'][i,0] > swl)[0][0]
            k = np.argwhere(data['vol'][i,1] < swl)[0][0]
            data['isw'][i] = data['cur'][i,0,j]
            print('Switching current %2.e at Vg = %.2f' % (data['isw'][i],ggv[i]))
            data['irt'][i] = data['cur'][i,1,k]
            print('Retrapping current %2.e at Vg = %.2f' % (data['irt'][i],ggv[i]))
        except IndexError:
            data['isw'][i] = 0
            print('No switching to normal state at Vg = %.2f (criteria %.2e V)' % (ggv[i],swl))
            data['irt'][i] = 0
            print('No retrapping to sc state at Vg = %.2f (criteria %.2e V)' % (ggv[i],swl))
    # Save to file
def plms(dev,opt,plp):
    """Switching current distribution measurement program."""

    # Importing modules
    from instruments import wfsg,mult,dcpw,rbrg         # import instrument classes
    from analysis import vol2ppr                        # import pulse counter funtion
    import numpy as np                                  # import numpy for array operations
    import PyDAQmx as nidaq                             # Import DAQ module
    import time                                         # For sleep()
    import matplotlib.pyplot as plt
    plt.ion()
    #from PyQt4 import QtGui
    from pyqtgraph.Qt import QtGui, QtCore
    import pyqtgraph as pg

    if 'wgm' in dev:                                    # set main waveform generator
        wgm = wfsg(dev['wgm'])
    else:
        raise ValueError('No main wfsg')

    if 'wgg' in dev:                                    # set gate generator
        wgg = wfsg(dev['wgg'])
    elif 'psg' in dev:
        wgg = dcpw(dev['psg'])
    else:
        raise ValueError('No gate generator')

    if 'brg' in dev:                                    # set resistance bridge
        brg = rbrg(dev['brg'][0],dev['brg'][1])
    else:
        raise ValueError('No resistance bridge')

    if 'nic' in dev:                                    # set NI DAQ card
        t = nidaq.Task()
        t.CreateAIVoltageChan(                          # Creating analog input voltage channel
                                "Dev1/ai0",             # Device name / recorded channel
                                None,
                                nidaq.DAQmx_Val_Diff,
                                0,                      # MIN amplitude
                                0.1,                    # MAX amplitude
                                nidaq.DAQmx_Val_Volts,  # Units
                                None
                                )

    # General settings
    bre = opt['bre']                                    # Bias resistance
    gnv = opt['gnv']                                    # gain for voltage amplifier
    gnc = opt['gnc']                                    # gain for current amplifier
    ggp = opt['ggp']                                    # gate generator pause
    swl = opt['swl']                                    # switching level for the jump to normal state
    gss = opt['gss']                                    # gate sweep speed
    ivb = opt['ivb']                                    # bias type
    gcs = opt['gcs']                                    # Gate on/off

    # Pulse settings
    frq_gen = plp['frq_gen']                            # Pulse main frequency = 1/(2*t)
    frq_bum = plp['frq_bum']                            # Pulse burst frequency = 1/T
    tfs = plp['tfs']                                    # Plot recorded pulses on/off
    gtc = plp['gtc']                                    # wgm trigger on/off
    amp = plp['amp']                                    # Current value of 1st pulse
    cps = plp['cps']                                    # Current amplitude step
    pst = plp['pst']                                    # Pulse statistics value
    ggv = plp['ggv']                                    # Gate voltage value

    # Setting up data
    data = {
            'cur': np.zeros(1),                         # Dynamic list of current values
            'spr': np.zeros(0)                          # Dynamic list of switching probabilities
            }
    data['cur'][0] = amp                         # Initialize 1st value
    # Journal

    # Save file

    # Figures
    app = QtGui.QApplication([])
    win = pg.GraphicsWindow(title="Scatter Plot Symbols")
    win.resize(1000,600)
    pg.setConfigOptions(antialias=True)
    plot = win.addPlot(title="Plotting with symbols")
    plot.addLegend()
    # Set up pulse & NI card
    wgm.set_pul_par(frq_gen,50,frq_bum,1,0)             # Setting wgm to single pulse mode
    prt = 1.2*(1/frq_bum)                               # Recording time for 1 pulse + interval
    t.CfgSampClkTiming(                                 # Setting up task for NI DAQ card
                        "",
                        int(1e6),                       # Samples per second
                        nidaq.DAQmx_Val_Rising,
                        nidaq.DAQmx_Val_FiniteSamps,
                        int(pst*prt*1e6)                # Measurement time x sample rate
                        )
    pvr = np.zeros((int(pst*prt*1e6),), dtype=np.float64)# Voltage responce array
    read = nidaq.int32()                                # Data type read
    # Read temperature
    # Journal
    # Figures

    # Start measurements
    try:
        wgg.set_gate(ggv,gss,gcs,ggp,ivb)               # Set gate voltage

        # Pulse cycle
        u = 0
        i = 0
        while u < 8:
            wgm.set_pul_amp(data['cur'][i]*bre)
            time.sleep(0.05)

            # Recoding the pulses
            t.StartTask()                                   # Start NI card recording
            # Sending triggers
            for j in range(pst):
                wgm.trg(gtc)
                time.sleep(1/frq_bum)
            t.ReadAnalogF64(                                # Read data from the NI card
                            int(pst*prt*1e6),               # Amount of data
                            25,                             # Delay if data is not ready
                            nidaq.DAQmx_Val_GroupByChannel,
                            pvr, len(pvr),                  # where to read and how much
                            nidaq.byref(read),              # buffer ?
                            None
                            )
            t.StopTask()                                    # Stop NI card recording
            # Analysis of recording
            plt.plot(pvr)
            ppr = vol2ppr(pvr,swl)/pst                      # Probability from recording analysis
            # Correcting data
            if ppr > 1:
                data['spr'] = np.append(data['spr'],1)
            else:
                data['spr'] = np.append(data['spr'],ppr)
            print('Probability %.3f at current %.3e' % (ppr,data['cur'][i]))
            if (ppr - 1) < 0.05:
                u = u + 1
            else:
                u = 0
            if u < 8:
                data['cur'] = np.append(data['cur'],data['cur'][i] + cps)
                i = i + 1
            else:
                break
        lep = len(data['cur'])                              # Length of recorded SCD
        mpi = round(lep/2)                                  # ~ middle point number
        wgm.set_pul_amp(data['cur'][mpi]*bre)               # set pulse amplitude
        time.sleep(0.05)
        t.StartTask()
        for j in range(pst):
            wgm.trg(gtc)
            time.sleep(1/frq_bum)
        t.ReadAnalogF64(                                # Read data from the NI card
                        int(pst*prt*1e6),               # Amount of data
                        15,                             # Delay if data is not ready
                        nidaq.DAQmx_Val_GroupByChannel,
                        pvr, len(pvr),                  # where to read and how much
                        nidaq.byref(read),              # buffer ?
                        None
                        )
        t.StopTask()

        #plt.plot(pvr)
        #plt.show()
        # Analysis of recording
        ppr = vol2ppr(pvr,swl)/pst
        print('Probability %.3f at current %.3e' % (ppr,data['cur'][mpi]))
        plot.plot(pvr, pen=(0,0,200), symbolBrush=(0,0,200), symbolPen='w', symbol='o', symbolSize=14, name="symbol='o'")
        #plt.plot(data['cur'],data['spr'])
        #plt.show()
        import sys
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()
    except:
        raise
    finally:
        # Closing everything
        wgg.set_gate(0,gss,gcs,0,ivb)
        wgm.set_bum(0)
        wgm.set_out(0)
        wgm.close()
        wgg.close()

    # Final temperature measurement
    # Save file
    # Final plot
