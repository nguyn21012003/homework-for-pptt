import numpy as np
import time
from multiprocessing import Pool

# Constants
pi = np.pi
hbar = 658.5  # mev fs
e = np.e
i = complex(0, 1)

# Parameters
a0 = 125  # Ban kinh Borh
khi0 = 0.0001  # Chi_0
T20 = 210  # Thoi gian hoi phuc
ga = 6.5 * 1e-20  # Lien quan mat do
deltat = 25  # Thoi gian gi do
emax = 300  # mev
Delta0 = 30  # mev Deturning
Er = 4.2  # mev
N = 100 # So phuong trinh

# Some predefined constants
deltae = emax / N
prefac = deltae * np.sqrt(deltae) * (1e24) / (pi**2 * Er**(3/2) * a0**3)
fpn = np.zeros((2, N + 1), np.complex_)

# Initial conditions
t0 = -3 * deltat
tm = 1000
h = 2
M = int((tm - t0) / h)

# Define essential functions
def g(n, l):
    domi = np.sqrt(n) - np.sqrt(l)
    check = domi == 0
    domi[check] = 1.0
    pref = 1.0 / np.sqrt(n * deltae)
    fac = np.abs((np.sqrt(n) + np.sqrt(l)) / domi)
    g = pref * np.log(fac)
    g[check] = 0.0
    return g

def OM(t, gp):
    pref1 = np.sqrt(pi) / (2. * deltat) * khi0 * np.exp(-t * t / (deltat * deltat))
    pref2 = np.sqrt(Er) / (hbar * pi) * deltae * gp
    return pref1 + pref2

def E(gf):
    return (np.sqrt(Er) / pi) * deltae * gf

# Optimized fp function with unrolled loops
def fp(t, fpn, T2):
    arrjj = np.arange(1, N + 1)
    fp = np.zeros((2, N + 1), np.complex_)
    
    for n in range(1, N + 1, 2):  # Loop incremented by 2
        arrj = np.delete(arrjj, n - 1)
        gf = np.sum(g(n, arrj) * 2 * fpn[0, arrj])
        gp = np.sum(g(n, arrj) * fpn[1, arrj])

        OMM = OM(t, gp)
        fp[0, n] = np.imag(-2. * (OMM * np.conjugate(fpn[1, n])))
        fp[1, n] = (-(i / hbar) * (n * deltae - Delta0 - E(gf)) - 1. / T2) * fpn[1, n] + i * (1. - 2. * fpn[0, n]) * OMM

        # Unroll next iteration
        if n + 1 <= N:
            arrj = np.delete(arrjj, n)
            gf = np.sum(g(n + 1, arrj) * 2 * fpn[0, arrj])
            gp = np.sum(g(n + 1, arrj) * fpn[1, arrj])

            OMM = OM(t, gp)
            fp[0, n + 1] = np.imag(-2. * (OMM * np.conjugate(fpn[1, n + 1])))
            fp[1, n + 1] = (-(i / hbar) * ((n + 1) * deltae - Delta0 - E(gf)) - 1. / T2) * fpn[1, n + 1] + i * (1. - 2. * fpn[0, n + 1]) * OMM

    return fp

# Song song hóa hàm rk4
def rk4_parallel(args):
    t, fpn, T2, start, end = args
    Nt_part = np.zeros(end - start, np.float_)
    Pt_part = np.zeros(end - start, np.complex_)
    tt_part = t + h * np.arange(start, end)
    
    for j in range(start, end):
        k1 = h * fp(t, fpn, T2)
        t_new = t + h / 2.
        k2 = h * fp(t_new, fpn + k1 / 2., T2)
        k3 = h * fp(t_new, fpn + k2 / 2., T2)
        k4 = h * fp(t_new + h / 2, fpn + k3, T2)
        
        fpn += (k1 + 2. * k2 + 2. * k3 + k4) / 6.
        
        Nt_part[j - start] = np.real(prefac * np.sum(fpn[0, :] * np.sqrt(np.arange(N + 1))))
        Pt_part[j - start] = prefac * np.sum(fpn[1, :] * np.sqrt(np.arange(N + 1)))
        t += h
        
        T2 = 1 / (1 / T20 + ga * Nt_part[j - start])
        
    return tt_part, Nt_part, Pt_part

# Hàm rk4 chính, sử dụng multiprocessing để gọi rk4_parallel song song
def rk4(t, fpn, T2, num_core):
    # Tính kích thước mỗi phần công việc (chunk) cho mỗi lõi
    chunk_size = M // num_core 
    chunks = []
    
    # Phân chia công việc thành các phần dựa vào num_core
    for core in range(num_core):
        start = core * chunk_size + 1
        end = (core + 1) * chunk_size + 1 if core < num_core - 1 else M + 1
        chunks.append((t, fpn.copy(), T2, start, end))
    
    # Sử dụng multiprocessing Pool để thực thi song song
    with Pool(num_core) as pool:
        results = pool.map(rk4_parallel, chunks)
    
    # Ghép kết quả từ các lõi
    tt = np.concatenate([result[0] for result in results])
    Nt = np.concatenate([result[1] for result in results])
    Pt = np.concatenate([result[2] for result in results])
    
    return tt, Nt, Pt

def Fourier(tt, Nt, Pt):
    step = 1000
    dt = h
    a, b = -100, 100
    deltapi = (b - a) / step
    earr = np.linspace(a, b, step)
    Pp = dt * np.array([np.sum(Pt * e**(i * (a + n * deltapi) * tt / hbar)) for n in range(step)])
    Ep = dt * np.array([np.sum(E(tt) * e**(i * (a + n * deltapi) * tt / hbar)) for n in range(step)])
    res1 = np.imag(Pp / Ep)
    # print(earr.shape,res1.shape)
    # exit()
    writetofileFT(earr, res1, np.real(Pp), np.imag(Pp), np.real(Ep), np.imag(Ep))
# Main function
def main():
    t = t0
    fpn[0, :] = np.zeros(N + 1, np.float_)
    fpn[1, :] = np.zeros(N + 1, np.complex_)
    set_cores = 7
    start_time = time.time()
    tt, Nt, Pt = rk4(t, fpn, T20,set_cores)
    
    writetofile(tt, Nt, Pt)
    Fourier(tt, Nt, Pt)    

    print(f"Process time = {time.time() - start_time:.2f} seconds")

# Write to file functions
def writetofileFT(t, res1, real_Pp, imag_Pp, real_Ep, imag_Ep):
    with open('./dataFT.txt', 'w') as c:
        for i in range(len(t)):
            c.write(f"{t[i]:e}\t{res1[i]:e}\t{real_Pp[i]:e}\t{imag_Pp[i]:e}\t{real_Ep[i]:e}\t{imag_Ep[i]:e}\n")
    print('Data written to dataFT.txt')

def writetofile(t, a, b):
    with open('./data.txt', 'w') as c:
        for i in range(1, len(a)):
            c.write(f"{t[i]:e}\t{a[i]:e}\t{np.abs(b[i]):e}\n")
    print('Data written to data.txt')

if __name__ == '__main__':
    main()
