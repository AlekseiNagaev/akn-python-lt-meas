import sys
import traceback
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rc('font', size=28)
from scipy.constants import pi, hbar, h
from math import exp, log, pow, sqrt, acos, cosh, tanh, sin
from scipy.optimize import curve_fit, fmin, fsolve, brentq
from scipy.integrate import quad
import nuf,conv,tef
try:
	str = '100'#input('File name?\n')
	str = "PL@" + str + "mK.mat"
	cur, pr, nT  = tef.load_file(str)
	T = conv.KeV(nT)
	x = cur
	print("Temperature is %.2f K / %.2e eV" % (nT,T))
	print("Critical temperature of the junction %.2f K / %.2e eV" % (nuf.nTc,nuf.Tc))
	print("Junction capacitance %.2e F" % nuf.C)
	print("Number of superconducting channels num/exp %d %e" % (nuf.Nch,nuf.eNch))
	#nuf.C = 1e-11
	print("Transmission coefficient %.3f" % nuf.tau)
	print("Normal state resistance %.2f Ohm" % nuf.Rn)
	print("Normal state resistivity %.2e Ohm*m" % nuf.rn)
	print("Charging energy %.2e J / %.2e eV" % (nuf.Ec0,nuf.Ec))
	print("Zero-temperature gap value %.2e eV" % nuf.d0 )
	print("Expected gap value is %.2e eV" % nuf.delta(T))
	phi = np.linspace(-4*pi, 4*pi, 1000)
	Ij = nuf.nIj(phi,T)
	U0 = nuf.nU(phi,0,T)
	Ic, phi0 = nuf.nIc(T)
	Ub0, p10, p20 = nuf.nUb(0,T,phi0)
	print(r'Critical current %.2e A Critical phase %.2f$\pi$' % (Ic,phi0/pi))
	wp = nuf.w0(phi0,T)
	print("Plasma frequency at critical phase %.2e Hz" % wp)
	nTp = conv.JK(conv.HzJ(wp))
	Tp = conv.JeV(conv.HzJ(wp))
	print("Approximate crossover temperature %.2f K / %.2e eV" % (nTp,Tp))
	Ib = np.linspace(0.8*Ic, 0.99*Ic, 1000)
	l = len(phi)
	j = len(Ib)
	Ub = np.zeros(j)
	#Gta = np.zeros(j)
	#Gqt = np.zeros(j)
	Pta = np.zeros(j)
	Pqt = np.zeros(j)
	for i in range(j):
		Ub[i], p1, p2 = nuf.nUb(Ib[i],T,phi0)
		if T>Tp:
			#Gta[i] = nuf.nGta(Ib[i],T)
			Pta[i] = nuf.nPta(Ib[i],T)
		else:
			#Gqt[i] = nuf.nGqt(Ib[i],
			Pqt[i] = nuf.nPqt(Ib[i],T)

	fig = plt.figure(figsize=(10, 9))
	t = 0
	if t:
		ax1 = fig.add_subplot(221)
		ax1.plot(phi/pi, Ij, 'bo', label='Ij')
		ax1.axvline(p10/pi,color='r')
		ax1.axvline(p20/pi,color='purple')
		ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g$\pi$'))
		ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0))
		ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

		ax2 = fig.add_subplot(223)
		ax2.plot(phi/pi, U0, 'bo', label='U0')
		ax2.axvline(p10/pi,color='r')
		ax2.axvline(p20/pi,color='purple')
		ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g$\pi$'))
		ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0))
		ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

		ax3 = fig.add_subplot(122)
	else:
		ax3 = plt.axes()

	ax3.grid(1)
	if T>Tp:
		#ax1.plot(Ib, Gta, 'bo', label='Gta')
		ax3.plot(Ib, Pta, 'bo', label='Pta')
	else:
		#ax1.plot(Ib, Gqt, 'bo', label='Gqt')
		ax3.plot(Ib, Pqt, 'bo', label='Pqt')
	#ax2 = fig.add_subplot(212)
	ax3.plot(cur,pr,"ro",label="P_exp")
	ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	ax3.legend()
	plt.xlabel('I, A', fontsize=32)
	plt.ylabel('P', fontsize=32)
	plt.savefig('test.png', bbox_inches='tight')
	#plt.show()
except:
	traceback.print_exc()
	print('Error\n',sys.exc_info()[0])
	#input("Press Enter to continue...")
	#raise
#input("Press Enter to continue...")
