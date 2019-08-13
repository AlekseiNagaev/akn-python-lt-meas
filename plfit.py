import sys
import traceback
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', size=28)
import scipy.constants as const
from math import exp, log, pow, sqrt, acos
from scipy.io import loadmat
from scipy.optimize import curve_fit
C = 135e-15 # Capacitance of the junction
t = 125e-6#4e3 # Length of the current pulses
Q = 5 # Quality factor of the junction
m1 = 1 # Reduction factor for RCSJ TA model func1
m2 = 0.7 # Reduction factor for cubic TA model func2
m3 = 1 # Reduction factor for RCSJ MQT model func3
m4 = 0.7 # Reduction factor for cubic MQT model func4
m5 = 1 # Reduction factor for cubic PD model func5
try:
	def wp0(b):
		return sqrt(2*b*const.e/(const.hbar*C))
		
	def wp(x,b):
		return sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(1-np.power(np.divide(x,b),2),0.25)
	
	def func1( x, a, b): # RCSJ TA model
		wp = sqrt(2*const.e*b/(const.hbar*C))*np.power(1-np.power(np.divide(x,b),2),0.25)
		dU = m1*2*(const.hbar/(2*const.e))*b*(np.sqrt(1 - np.power(np.divide(x,b),2)) - np.divide(x,b)*np.arccos(np.divide(x,b)))
		at = sqrt(1 + 1/(4*Q^2)) - 1/(2*Q) 
		G = at*wp/(2*const.pi)*np.exp(-dU/(const.k*a))
		return 1 - np.exp(-G*t)

	def func2( x, a, b): # Cubic TA model
		wp = sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(2*(1-np.divide(x,b)),0.25)
		dU = m2*0.66*(const.hbar/(2*const.e))*b*sqrt(8)*np.power((1-np.divide(x,b)),1.25)
		at = sqrt(1 + 1/(4*Q**2)) - 1/(2*Q)
		G = at*wp/(2*const.pi)*np.exp(-dU/(const.k*a))
		return 1 - np.exp(-G*t)

	def func3( x, a): # RCSJ MQT model
		wp = sqrt(2*const.e/(const.hbar*C))*sqrt(a)*np.power(1-np.power(np.divide(x,a),2),0.25)
		dU = m3*2*(const.hbar/(2*const.e))*a*(np.sqrt(1 - np.power(np.divide(x,a),2)) - np.divide(x,a)*np.arccos(np.divide(x,a)))
		aq = np.sqrt(120*const.pi*7.2*dU/(const.hbar*wp))
		G = aq*wp/(2*const.pi)*np.exp(-7.2*dU/(const.hbar*wp)*(1+0.87/Q))
		return 1 - np.exp(-G*t)

	def func4( x, a): # Cubic MQT model
		wp = sqrt(2*const.e/(const.hbar*C))*sqrt(a)*np.power(2*(1-np.divide(x,a)),0.25)
		dU = m4*0.66*(const.hbar/(2*const.e))*a*sqrt(8)*np.power((1-np.divide(x,a)),1.5)
		G = 12*wp*sqrt(6*const.pi)/(2*const.pi)*np.sqrt(np.divide(dU,wp))/sqrt(const.hbar)*np.exp(-7.2*dU*(1+0.87/Q)/(const.hbar*wp))
		return 1 - np.exp(-G*t)
	
	def func5 (x, a, b, c): # Cubic PD model
		Ej = const.hbar*b/(2*const.e) 
		wp = sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(2*(1-np.divide(x,b)),0.25)
		dU = m5*0.66*Ej*sqrt(8)*np.power((1-np.divide(x,b)),1.25)
		dUrt = Ej*((Q**2)/2)*(((x - c)**2)/(b**2))
		Grt = wp*((x - c)/b)*sqrt(Ej/(2*const.pi*const.k*a))*np.exp(-dUrt/(const.k*a))
		Prt = 1 - np.exp(-Grt*t)
		at = sqrt(1 + 1/(4*Q**2)) - 1/(2*Q)
		Gta = at*wp/(2*const.pi)*np.exp(-dU/(const.k*a))
		Gs = Gta*(1-np.divide(1,Prt))*np.log(1-Prt)
		return 1 - np.exp(-Gs*t)
	
	def prt(x, a, b, c):
		Ej = const.hbar*b/(2*const.e) 
		wp = sqrt(2*const.e/(const.hbar*C))*sqrt(b)*np.power(2*(1-np.divide(x,b)),0.25)
		dUrt = Ej*((Q**2)/2)*(((x - c)**2)/(b**2))
		Grt = wp*((x - c)/b)*sqrt(Ej/(2*const.pi*const.k*a))*np.exp(-dUrt/(const.k*a))	
		return 1 - np.exp(-Grt*t)
		
	def load_file( str ):
		mdata = loadmat(str)
		data = mdata['data']
		if 'pul_i' in data.dtype.names:
			cur = mdata['data']['pul_i'][0,0][0].tolist()
		elif 'cur' in data.dtype.names:
			cur = mdata['data']['cur'][0,0][0].tolist()
		else:
			print("Exception #1.\n Peace out!")
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

	str = '600'#input('File name?\n')
	str = "PL@" + str + "mK.mat"
	cur, pr, T  = load_file(str)
	x = cur
	y = pr
	plt.figure(figsize=(10,9))
	ed = 1
	if ed:
		plt.plot(x, y, 'bo', label='data')
	ta = 1
	qt = 1
	if (T > 0.119) and ta: #crossover temperature
		po1, pc2 = curve_fit(func1, x, y, bounds=([T, 19.5e-7], [T + 1e-3, 22e-7]))
		po2, pc2 = curve_fit(func2, x, y, bounds=([T, 5.42e-7], [T + 1e-3, 5.82e-7]))
		a = po1[0]
		b = po1[1]
		a2 = po2[0]
		b2 = po2[1]
		print('fit: a=%1.2e, b=%1.2e' % (a, b))
		wp1 = wp0(b)
		print('TA $w_{p0}$ = %.2e' % wp1)
		print('fit2: a=%1.2e, b=%1.2e' % (a2, b2))
		wp2 = wp0(b2)
		print('TA3 $w_{p0}$ = %.2e' % wp2)
		#plt.plot(x, func1(x, *po1), 'r--', label='TA fit: $T_m = %.2f$, $I_m$=%1.2e' % (a, b))
		plt.plot(x, func2(x, *po2), 'g--', label='TA3 fit: $T_m = %.2f$, $I_m$=%1.2e' % (a2, b2))
	elif qt:
		po3, pc3 = curve_fit(func3, x, y, bounds=([17e-7], [21e-7]))
		po4, pc4 = curve_fit(func4, x, y, bounds=([4.73e-7], [6.11e-7]))
		a3 = po3[0]
		a4 = po4[0]
		#plt.plot(x, func3(x, *po3), 'm--', label='MQT fit: $I_a$=%1.2e' % a3)
		print('wp0 = %1.2e' % wp0(a3))
		Tc3 = const.hbar*wp0(a3)/(7.2*const.k)
		print('Tc = %1.3f' % Tc3)
		plt.plot(x, func4(x, *po4), 'g--', label='MQT3 fit: $I_a$=%1.2e' % a4)
		print('wp0 = %1.2e' % wp0(a4))
		Tc4 = const.hbar*wp0(a4)/(7.2*const.k)
		print('Tc = %1.3f' % Tc4)
	m = 0
	if m:
		am =  0.3
		bm = 18e-7
		plt.plot(x, func1(x, am , bm), 'g-o', label='TA mfit: $T_m = %.2f$, $I_m$=%1.2e' % (am, bm))
		plt.plot(x, func3(x, bm), 'g-o', label='MQT mfit: $I_m$=%1.2e' % bm)
	cu = 0
	if cu:
		Im = 18e-7 
		x2 = np.linspace(2e-6, 1e-4, num=50000)
		plt.plot(x2, func3(x2, Im), 'g-o', label='man_fit: $I_m$=%1.2e' % Im)	
	rt = 0
	if rt:
		a = T
		b = 4.15e-7
		c = 1.3e-7
		plt.plot(x, prt(x, a, b, c), 'y--', x, func5(x, a, b, c),'r--')
	pd = 0
	if pd:
		po3, pc3 = curve_fit(func5, x, y, bounds=([T, 4.2e-7, 2.1e-7], [T + 1e-3, 6.5e-7, 2.6e-7]))
		a3 = po3[0]
		b3 = po3[1]
		c3 = po3[2]
		plt.plot(x, func5(x, *po3), 'm--', label='PD3 fit: $T_m = %.2f$, $I_m$=%1.2e, $I_r$=%1.2e' % (a3, b3, c3))
	plt.xlabel('I, A', fontsize=32)
	plt.ylabel('P', fontsize=32)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	plt.tick_params(axis='both', which='major', labelsize=32)
	plt.legend(fontsize=20)
	plt.show()
except:
	traceback.print_exc()
	#print('Error\n',sys.exc_info()[0])
	#raise
input("Press Enter to continue...")
