from scipy.io import loadmat
import numpy as np
def load_file(str):
  mdata = loadmat(str)
  data = mdata['data']
  if 'pul_i' in data.dtype.names:
    cur = mdata['data']['pul_i'][0,0][0].tolist()
  elif 'cur' in data.dtype.names:
    cur = mdata['data']['cur'][0,0][0].tolist()
  else:
    print("Exception #1.\n Wrong data names in file!")
    sys.exit()
  if len(cur) == 1:
    cur = list(map(list, zip(*cur)))
  pr = mdata['data']['pr'][0,0].tolist()
  if len(pr) != 1:
    pr = list(map(list, zip(*pr)))
    pr = pr[0]
  T = mdata['T1'][0,0]
  cur = np.asarray(cur)
  pr = np.asarray(pr)
  T = np.asarray(T)
  return cur, pr, T
