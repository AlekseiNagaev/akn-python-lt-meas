def vol2ppr(vol,swl):
    import numpy as np
    arr = np.where(vol>swl)
    arr = np.diff(arr)
    arr = np.equal(arr,1)
    sfl = False
    cou = 0
    for i in range(len(arr[0])):
        if arr[0][i]:
            if not sfl:
                cou = cou + 1
                sfl = True
        else:
            sfl = False
    return cou
