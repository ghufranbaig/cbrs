import math
import sys
d = 1

def dist(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.hypot(x,y)

def angle(a,b):
	x = a[0] - b[0]
	y = a[1] - b[1]
	return math.atan2(y,x)

def dBmtoW(dbm):
	return (10**((dbm-30)/10.0))

def WtodBm(W):
	return (10*math.log10(W) + 30)

def RxPowertoDistance(P_r,P_t):
	n = 3.5
	L_0 = 25.686
	d_0 = 1
	return( d_0 * 10**( math.log10(P_t/P_r)/n - L_0/(n*10.0) ) )

def RxPwrtoTxPower (P_r,d):
	n = 3.5
	L_0 = 25.686
	d_0 = 1
	return ( P_r*(10**((L_0 + 10*n*math.log10(d/d_0))/10.0)) )

def getOverHearingRange(d):
	#d = dist(enb,ue)
	TxP = RxPwrtoTxPower (dBmtoW(-90.0),d)
	m_range = RxPowertoDistance(dBmtoW(-100.0),TxP)
	return m_range


if len(sys.argv) > 1:
   d = int(sys.argv[1])

#d *= 2.78

print WtodBm(RxPwrtoTxPower (dBmtoW(-95),d))
print RxPowertoDistance(dBmtoW(-99),dBmtoW(30))
P_t = 7 # 29 dbm

n = 3.5
L_0 = 25.686
d_0 = 1
P_r = P_t/(10**((L_0 + 10*n*math.log10(d/d_0))/10))

D = d_0 * 10**( math.log10(P_t/P_r)/n - L_0/(n*10) )

P_T = P_r*(10**((L_0 + 10*n*math.log10(d/d_0))/10))


P = P_r
dbm = 10*math.log10(P/0.001)

#print(D)
#print(P_T)

#print(P)
#print(dbm)
