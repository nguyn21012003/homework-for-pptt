#-------------------------------------------------------------------------------------------------------------------------------------#
import numpy as np
import time
#-------------------------------------------------------------------------------------------------------------------------------------#

#Constants:
#-------------------------------------------------------------------------------------------------------------------------------------#
pi      = np.pi
hbar    = 658.5 #mev fs
e       = np.e
i       = complex(0,1)

#Parameters:
khi0    = 0.1 
T2      = 100 #fs
deltat  = 100 #25 #fs
emax    = 300 #mev
Delta0  = 100#30 #mev
Er      = 4.2 #mev
N       = 50
#Some predefined constants
deltae  = emax/N
pref    = deltae*np.sqrt(deltae)
fpn     = np.zeros((2,N+1),np.complex_)
#-------------------------------------------------------------------------------------------------------------------------------------#


#Initial conditions:
#-------------------------------------------------------------------------------------------------------------------------------------#
t0      = -3*deltat

tm      = 500
h       = 2
M       = int((tm - t0)/h)
#-------------------------------------------------------------------------------------------------------------------------------------#


# Function that open and write to file
#-------------------------------------------------------------------------------------------------------------------------------------#
def writetofile(t,a,b): 
    c = open('./data.txt','w+')
    for i in range(1,M+1):
        c.write("%f \t %f \t %f \n" % (t[i],a[i],b[i])) 

    c.close()
    print('All most done! Congratulation !!?')
#Define essential functions:
#-------------------------------------------------------------------------------------------------------------------------------------#
def g(n,l):
    pref    = 1./np.sqrt(n*deltae)
    fac     = np.abs((np.sqrt(n) + np.sqrt(l))/(np.sqrt(n) - np.sqrt(l)))
    g       = pref*np.log(fac)
    return g
#-------------------------------------------------------------------------------------------------------------------------------------#   
def OM(t,gp):
    pref1   = np.sqrt(pi)/(2.*deltat)*khi0*np.exp(-t*t/(deltat*deltat))
    pref2   = np.sqrt(Er)/(hbar*pi)*deltae*gp
    om      = pref1 + pref2
    return om
#-------------------------------------------------------------------------------------------------------------------------------------#
def E(gf):

    E = ((np.sqrt(Er))/pi)*deltae*gf
    return E
#-------------------------------------------------------------------------------------------------------------------------------------#
def fp(t,fpn):
    arrjj = np.arange(1,N+1)
    fp = np.zeros((2,N+1),np.complex_)
    for n in range(1,N+1,2): 
    	#Inorging one loop here by replacing it by np.sum and np.arange
        arrj = np.delete(arrjj,n-1)
        gf = np.sum(g(n,arrj)*2*fpn[0,arrj])
        gp = np.sum(g(n,arrj)*fpn[1,arrj])
        #----------------------------------#
        
        OMM = OM(t,gp)
        fp[0,n] = np.imag(-2.*(OMM*np.conjugate(fpn[1,n])))
        fp[1,n] = (-(i/hbar)*(n*deltae - Delta0 -E(gf)) - 1./T2)*fpn[1,n] + i*(1.-2.*fpn[0,n])*OMM
        
        #Loop 2
        arrj = np.delete(arrjj,n)
        gf = np.sum(g(n+1,arrj)*2*fpn[0,arrj])
        gp = np.sum(g(n+1,arrj)*fpn[1,arrj])
        
        OMM = OM(t,gp)
        fp[0,n+1] = np.imag(-2.*(OMM*np.conjugate(fpn[1,n+1])))
        fp[1,n+1] = (-(i/hbar)*(n*deltae - Delta0 -E(gf)) - 1./T2)*fpn[1,n+1] + i*(1.-2.*fpn[0,n+1])*OMM
        #Faster loop by increasing step by 2, athough by 1 in the common sense !!!
        
        
    return fp  
#-------------------------------------------------------------------------------------------------------------------------------------#	
#Main rk4 :
def rk4(t,fpn):
    Nt = np.zeros(M+1,np.float_)
    Pt = np.zeros(M+1,np.float_)
    #Time array
    tt    = t + h*np.arange(0,M+1)
    for j in np.arange(1,M+1):
        k1=h*fp(t,fpn)
        t_new = t + h/2.
        k2=h*fp(t_new,fpn+k1/2.)
        k3=h*fp(t_new,fpn+k2/2.)
        k4=h*fp(t_new+h/2,fpn+k3)
    
        fpn = fpn+(k1 + 2.*k2 + 2.*k3 + k4)/6.

        Nt[j] = np.real(pref*np.sum(fpn[0,:]*np.sqrt(np.arange(N+1))))
        Pt[j] = np.abs(pref*np.sum(fpn[1,:]*np.sqrt(np.arange(N+1))))
        t = t + h
    writetofile(tt,Nt,Pt)
#-------------------------------------------------------------------------------------------------------------------------------------#

def main():
    t       = t0
    fpn[0,:]= np.zeros(N+1,np.float_)
    fpn[1,:]= np.zeros(N+1,np.complex_)
    begin_t = time.time()
    rk4(t,fpn)
    end = time.time()
    print("process time = %f (s) \n" % (end - begin_t))

if __name__ == '__main__':
    main()
