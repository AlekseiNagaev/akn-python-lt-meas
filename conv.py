from scipy.constants import e,k,hbar,h
def JeV(x):
	return x/e

def eVJ(x):
	return x*e

def eVK(x):
	return x/8.62e-5

def KeV(x):
	return x*8.62e-5

def KJ(x):
	return k*x

def JK(x):
	return x/k

def HzJ(x):
	return h*x

def JHz(x):
	return x/h
