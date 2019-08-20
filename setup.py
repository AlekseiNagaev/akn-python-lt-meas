def set_opt():

    # Import numpy for gate voltage
    import numpy as np

    # Device addresses
    dev = {
            'wgm': 12,
            'wgg': 13,
            'mu1': 10,
            'mu2': 11,
            'brg': [20,1],
            'lia': 14,
            'nic': 'Dev1'
            }
    # General options
    opt = {
            'bre': 1e6,
            'gnv': 1e2,
            'gnc': 1e5,
            'ggp': 0.05,
            'swl': 0.01,
            'gss': 0.05,
            'ivb': 'vol',
            'gcs': 1
            }
    # IV parameters
    ivp = {
            'ivs': 0.01,
            'pol': 0,
            'amp': 1,
            'stp': 0.05,
            'ggv': [0]#np.concatenate((np.arange(-10,10,0.5),10))
            }
    # Pulse parameters
    plp = {
            'frq_gen': 4e3,
            'frq_bum': 5,
            'tfs': 1,
            'gtc': 1,
            'amp': 1e-6,
            'cps': 1e-7,
            'pst': int(1e1),
            'ggv': 0
            }
    yield dev
    yield opt
    yield ivp
    yield plp
