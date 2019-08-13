import scipy.constants as const
from scipy.constants import pi, h, hbar,e
import conv
import numpy as np
from scipy.optimize import curve_fit, fmin, fsolve, brentq
from scipy.integrate import quad
from math import exp, log, pow, sqrt, acos, cosh, tanh, sin
from scipy.signal import argrelextrema

# Numerical calculations
t = 125e-6#4e3 # Length of the current pulses [s]

nTc = 1.32# [K]

Tc = conv.KeV(nTc)# [eV]

C = 344e-15 # [F]

Nch = 66.2 # 1 ch ~

tau = 0.9 # transmission coefficient

dp = 1e-4 #	delta for diff

Rq = h/e**2 # Resistance quantum 25812 [Ohm]

Rn = 100 # Normal resistance [Ohm]

w= 4e-6			# Width of graphene
h = 3.4e-10		# Thickness of graphene
l = 325e-9		# Length of graphene

rn = Rn*h*w/l

eNch = rn*pi*hbar/e**2

Ec0 = e**2/(2*C)# Charging energy in SI [J]

Ec = conv.JeV(Ec0) #value ~ 0.594e-6 [eV]

d0 = 1.764*Tc*80/100 # Gap value [eV]


def delta(T):# eV
	#print(T)
	return d0*np.sqrt(1 - np.power((T/Tc),3.03))						#eV
def nIj(phi,T):# A
	s0 = np.sqrt(1-tau*np.power(np.sin(phi/2),2))						#1
	nIj0 = delta(T)*s0/(2*T)											#1
	nIj0 = np.tanh(nIj0)												#1
	nIj0 = np.divide(nIj0,Rq*s0)										#1/Ohm
	nIj0 = Nch*tau*pi*delta(T)*np.multiply(nIj0,np.sin(phi))			#A
	#print("Josephson current %.1f A" % nIj0)
	return nIj0
def ndIj(phi,T): # A*s
	nd = (nIj(phi+dp,T)-nIj(phi,T))/dp
	#print("Current derivative %.1f" % nd)
	return nd
def nIc(T):# A
	def ndIj0(phi):
		return ndIj(phi,T)
	phi0 = fsolve(ndIj0,dp)
	return nIj(phi0,T), phi0
def w0(phi,T): #Hz
	w0 = np.sqrt(4*Ec*Rq*(nIj(phi+dp,T)-nIj(phi-dp,T))/(2*dp)) #eV
	#w0 = np.sqrt(Ec)
	w0 = w0*e/hbar
	#print("Plasma frequency at phase %.2f\u03C0 and temperature %.2f is %.2e Hz" % (phi/pi,conv.eVK(T),w0))
	#sqrt(dI/dphi)*sqrt(Ec/2*C)
	return w0
def nU(phi,Ib,T):# eV
	U0 = np.power(np.sin(phi/2),2)
	#print(U0)
	U0 = np.sqrt(1-tau*U0)
	#print(U0)
	U0 = delta(T)*U0/(2*T)
	#print(U0)
	U0 = np.cosh(U0)
	#print(U0)
	U0 = np.log(U0)
	#print(U0)
	U0 = -2*Nch*T*U0 - Ib*Rq*phi/(4*pi)	      #eV
	return U0

def nUb(Ib,T,phi0):# eV
	def nU0(phi,T):
		return nIj(phi,T)-Ib
	#phi1 = argrelextrema(x, np.less)
	#l = len(phi1)
	#l = l // 2
	#phi1 = phi1[l]
	#phi2 = argrelextrema(x, np.greater)
	phi1 = fsolve(nU0,dp,args=(T))
	phi2 = fsolve(nU0,pi,args=(T))
	nUb0 = nU(phi2,Ib,T) - nU(phi1,Ib,T)
	return nUb0, phi1, phi2

def nGta(Ib,T):
	phi0 = nIc(T)[1]
	#print("phi0 %.2f" % phi0)
	Ub, phi1, phi2 = nUb(Ib,T,phi0)
	#print("Ub %.2e eV" % Ub)
	#print("T %.2e eV" % T)
	nat = np.divide(Rq*Ec0,2*pi*Rn*hbar*w0(phi1,T))	#1
	#print("nat %.2e" % nat)
	nat = (np.sqrt(1 + np.power(nat,2)) - nat) 					#1
	#print("nat %.2e" % nat)
	Gta = np.exp(-Ub/T)									#1
	#print("Gta %.2e" % Gta)
	Gta = w0(phi1,T)*Gta/(2*pi)							#Hz
	#print("Gta %.2e" % Gta)
	Gta = nat*Gta												#Hz
	#print("Gta %.2e" % Gta)
	return Gta
def nGqt(Ib,T):
	phi0 = nIc(T)[1]
	Ub = nUb(Ib,T,phi0)[0]
	phi1 = nUb(Ib,T,phi0)[1]
	phi2 = nUb(Ib,T,phi0)[2]
	naq = np.divide(2.86*Rq*(conv.eVJ(Ec)*1e-6),2*pi*Rn*hbar*w0(phi1,T))
	#print("naq %.2e" % naq)
	naqt = 1 + 2.86*naq
	#print("naqt %.2e" % naqt)
	def dU(phi,Ib,T):
		return np.sqrt(nU(phi,Ib,T)-nU(phi1,Ib,T))
	dS = quad(dU,phi1,phi2,args=(Ib,T))
	#print("dS %e %e" % dS)
	S = (1 + 0.2777*naq)/sqrt(Ec)*dS[0]
	#print("S %.2e" % S)
	Gqt = np.exp(-S)
	#print("Gqt %.2e" % Gqt)
	Gqt = np.multiply(np.sqrt(120*pi*S),Gqt)
	#print("Gqt %.2e" % Gqt)
	Gqt = naqt*w0(phi1,T)*Gqt/(2*pi)
	#print("Gqt %.2e" % Gqt)
	return Gqt

def nPta(Ib,T):
	return 1-np.exp(-nGta(Ib,T)*t)
def nPqt(Ib,T):
	return 1-np.exp(-nGqt(Ib,T)*t)

#Fitting of experimental results
#C = 135e-15 # Capacitance of the junction
#Q = 5 # Quality factor of the junction
#m1 = 0.7 # Reduction factor for cubic TA model
#m2 = 0.7 # Reduction factor for cubic MQT model

#def wp0(b):
	#return sqrt(2*b*const.e/(const.hbar*C))

#def wp(x,b):
	#return sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(1-np.power(np.divide(x,b),2),0.25)

#def pta( x, a, b): # Cubic TA model
	#wp = sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(2*(1-np.divide(x,b)),0.25)
	#dU = m1*0.66*(const.hbar/(2*const.e))*b*sqrt(8)*np.power((1-np.divide(x,b)),1.25)
	#at = sqrt(1 + 1/(4*Q**2)) - 1/(2*Q)
	#Gta = at*wp/(2*const.pi)*np.exp(-dU/(const.k*a))
	#return 1 - np.exp(-Gta*t)

#def pqt( x, a): # Cubic MQT model
	#wp = sqrt(2*const.e/(const.hbar*C))*sqrt(a)*np.power(2*(1-np.divide(x,a)),0.25)
	#dU = m2*0.66*(const.hbar/(2*const.e))*a*sqrt(8)*np.power((1-np.divide(x,a)),1.5)
	#Gqt = 12*wp*sqrt(6*const.pi)/(2*const.pi)*np.sqrt(np.divide(dU,wp))/sqrt(const.hbar)*np.exp(-7.2*dU*(1+0.87/Q)/(const.hbar*wp))
	#return 1 - np.exp(-Gqt*t)
