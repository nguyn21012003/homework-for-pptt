#-------------------------------------------------------------------------------------------------------------------------------------#
import numpy as np
import time
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------------------------------------------------------------------#

#Constants:
#-------------------------------------------------------------------------------------------------------------------------------------#
pi      = np.pi
hbar    = 658.5 #mev fs
e       = np.e
i       = complex(0,1)

#Parameters:
a0      = 125           #Ban kinh Borh
khi0    = 0.0001            #Chi_0 
T20     = 210           #Thoi gian hoi phuc 
ga      = 6.5*1e-20     #Lien quan mat do 
deltat  = 25            #Thoi gian gi do 
emax    = 300           #mev
Delta0  = 30            #mev Deturning 
Er      = 4.2           #mev
N       = 100          # So phuong trinh 

#Some predefined constants
deltae     = emax/N
prefac      = deltae*np.sqrt(deltae)*(1e24)/(pi**2*Er**(3/2)*a0**3)
fpn         = np.zeros((2,N+1),np.complex_)
#-------------------------------------------------------------------------------------------------------------------------------------#


#Initial conditions:
#-------------------------------------------------------------------------------------------------------------------------------------#
t0      = -3*deltat

tm      = 1000
h       = 2
M       = int((tm - t0)/h)
print(M)
#-------------------------------------------------------------------------------------------------------------------------------------#


# Function that open and write to file
#-------------------------------------------------------------------------------------------------------------------------------------#
def writetofileFT(t,k,h,a,b,d): 
    c = open('./dataFT.txt','w+')
    for i in range(1,len(a)):
        c.write("%e \t %e \t %e \t %e \t %e \t %e \n" % (t[i],k[i],h[i],a[i],b[i],d[i]))

    c.close()
    print('All most done! Congratulation !!?')

def writetofile(t,a,b): 
    c = open('./data.txt','w+')
    for i in range(1,len(a)):
        c.write("%e \t %e \t %e  \n" % (t[i],a[i],np.abs(b[i])))

    c.close()
    print('All most done! Congratulation !!?')
#Define essential functions:
#-------------------------------------------------------------------------------------------------------------------------------------#
def g(n,l):
    domi    = np.sqrt(n) - np.sqrt(l)
    check   = np.where(domi==0.0)
    domi[check[0]] = 1.0
    pref    = 1./np.sqrt(n*deltae)
    fac     = np.abs((np.sqrt(n) + np.sqrt(l))/domi)
    g       = pref*np.log(fac)
    g[check[0]] = 0.0
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
def fp(t,fpn,T2):
    arrjj = np.arange(1,N+1)
    fp = np.zeros((2,N+1),np.complex_)
    gf = 0.0 
    gp = 0.0
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
        fp[1,n+1] = (-(i/hbar)*((n+1)*deltae - Delta0 -E(gf)) - 1./T2)*fpn[1,n+1] + i*(1.-2.*fpn[0,n+1])*OMM
        
        #Faster loop by increasing step by 2, athough by 1 in the common sense !!!

    return fp  
#-------------------------------------------------------------------------------------------------------------------------------------#	
#Main rk4 :

def rk4(t,fpn,T2):
    Nt = np.zeros(M+1,np.float_)
    Pt = np.zeros(M+1,np.complex_)
    #Time array
    tt    = t + h*np.arange(0,M+1)
    for j in np.arange(1,M+1):

        k1=h*fp(t,fpn,T2)
        t_new = t + h/2.
        k2=h*fp(t_new,fpn+k1/2.,T2)
        k3=h*fp(t_new,fpn+k2/2.,T2)
        k4=h*fp(t_new+h/2,fpn+k3,T2)
    
        fpn = fpn + (k1 + 2.*k2 + 2.*k3 + k4)/6. 

        Nt[j] = np.real(prefac*np.sum(fpn[0,:]*np.sqrt(np.arange(N+1))))
        Pt[j] = prefac*np.sum(fpn[1,:]*np.sqrt(np.arange(N+1)))
        t = t + h
        #T2 changes in time. 

        T2 = 1/(1/T20 + ga*Nt[j])


    return tt,Nt,Pt
#-------------------------------------------------------------------------------------------------------------------------------------#
def Et(t):
    res =  khi0*np.exp(-t*t/(deltat*deltat))
    return res

def Fourier(tt,Nt,Pt):
    step = 1000
    # Np = np.zeros((step),np.complex_)
    Pp = np.zeros((step),np.complex_)
    Ep = np.zeros((step),np.complex_)
    res1 = np.zeros((step))

    dt = h
    b = 100
    a = -100
    deltapi = (b-a)/step
    earr = np.linspace(a,b,step)
    # for n in range(1,step):
        # Ep[n] = Ew(earr[n],tt)
        # Pp[n] = Pw(earr[n],tt,Pt)
        # res1[n] = np.imag(Pp[n]/Ep[n])
    for n in np.arange(1,step):
        Pp[n] = dt*np.sum(Pt*e**(i*(a + n*deltapi)*tt/hbar))
        Ep[n] = dt*np.sum(Et(tt)*e**(i*(a + n*deltapi)*tt/hbar))
        res1[n] = np.imag(Pp[n]/Ep[n])
        # Np[n] = dt*np.sum(Nt)
        # Pp[n] = dt*np.sum(Pt)
        # Ep[n] = dt*np.sum(Et(a + deltapi*n))
        # res1[n] = np.imag(Pp[n]/Ep[n])
        # res2[n] = np.imag(Np[n]/Ep[n])

    writetofileFT(earr,res1,np.real(Pp),np.imag(Pp),np.real(Ep),np.imag(Ep))

def main():
    t       = t0
    fpn[0,:]= np.zeros(N+1,np.float_)
    fpn[1,:]= np.zeros(N+1,np.complex_)
    begin_t = time.time()
    tt,Nt,Pt = rk4(t,fpn,T20)
    writetofile(tt,Nt,Pt)
    Fourier(tt,Nt,Pt)
    end = time.time()
    print("process time = %f (s) \n" % (end - begin_t))

if __name__ == '__main__':
    main()
