import numpy as np
import matplotlib.pyplot as plt

class Ising:
    def make_all(self,L):
        res = Ising(L)
        res.generator_s(0)
        res.calc()
        res.plot_mag()
        res.plot_spec()
        return res
        

    def __init__(self,L,intv = 25,J = 1,B = 0):
        self.L = L
        self.N = L*L
        self.intv = intv
        self.T = [(5/self.intv)*(k+1) for k in range(self.intv)]
        self.J = J
        self.B = B
        self.cnt = 0
        self.ZU = np.zeros(intv,dtype=np.double)
        self.Z = np.zeros(intv,dtype=np.double)
        self.E2 = np.zeros(intv,dtype=np.double)
        self.E = np.zeros(intv,dtype=np.double)
        self.sc = np.ones(self.N,dtype=np.int0)

    def sweep(self): #Non-periodic
        self.cnt += 1
        res = 0
        for i in range(self.N):
            sum = 0

            nn = i + 1
            if(nn % self.L != 0):
                sum += self.sc[nn]
                # print(nn,sum)

            nn = i - 1
            if((nn+1) % self.L != 0):
                sum += self.sc[nn]
                # print(nn,sum)

            nn = i + self.L
            if(nn < self.N):
                sum += self.sc[nn]
                # print(nn,sum)

            nn = i - self.L
            if(nn >= 0):
                sum += self.sc[nn]
                # print(nn,sum)

            res += self.J*sum*self.sc[i]

        sigma = np.sum(self.sc)
        HH = -int(res/2) -self.B*sigma
        print(sigma,HH)
        for k in range(self.intv):
            T = (5/self.intv)*(k+1)
            ex = np.exp(-HH/T)
            self.ZU[k] += np.absolute(sigma)*ex
            self.Z[k] += ex
            self.E[k] += HH*ex
            self.E2[k] += HH*HH*ex
    
    def sweep2(self): # periodic
        self.cnt += 1
        res = 0
        # print(self.sc)
        for i in range(self.N):
            sum = 0

            nn = i + 1
            if(nn % self.L == 0): nn -= self.L
            sum += self.sc[nn]
                # print(nn,sum)

            nn = i + self.L
            if(nn >= self.N): nn -= self.N
            sum += self.sc[nn]
                # print(nn,sum)

            res += self.J*sum*self.sc[i]

        sigma = np.sum(self.sc)
        HH = -res -self.B*sigma
        # print(sigma,HH)
        for k in range(self.intv):
            T = (5/self.intv)*(k+1)
            ex = np.exp(-HH/T)
            self.ZU[k] += np.absolute(sigma)*ex
            self.Z[k] += ex
            self.E[k] += HH*ex
            self.E2[k] += HH*HH*ex

    def generator(self,i=0): # not using symmetry
        func = self.sweep2
        if(i == self.N): return
        if(i == 0): func()
        self.sc[i] *= -1
        func()
        self.generator(i+1)
        self.sc[i] *= -1
        self.generator(i+1)
        return
    
    def generator_s(self,i=1): # not using symmetry
        func = self.sweep2
        if(i == self.N): return
        if(i == 1): func()
        self.sc[i] *= -1
        func()
        self.generator(i+1)
        self.sc[i] *= -1
        self.generator(i+1)
        return
    
    def calc(self):
        self.mag = self.magnetization()
        self.spec = self.specific_heat()
        

    def magnetization(self): # not per spin
        res = np.zeros(self.intv)
        for i in range(self.intv):
            res[i] = self.ZU[i]/self.Z[i]
        return res
    
    def specific_heat(self):
        res = np.zeros(self.intv)
        for i in range(self.intv):
            res[i] = self.E2[i]/self.Z[i] - (self.E[i]/self.Z[i])**2
        return res

    def plot_mag(self, s= "save"):
        y = self.mag/max(self.mag)

        plt.style.use('seaborn-whitegrid')
        plt.ylim(0,1.1)
        plt.xlim(0,5)
        plt.plot(self.T,y)

        plt.ylabel('Abosolute magnetization per spin m')
        plt.xlabel('Temperature T')
        plt.show()
    
    def plot_spec(self,s ="save"):
        y = self.spec/max(self.spec)

        plt.style.use('seaborn-whitegrid')
        plt.ylim(0,1.1)
        plt.xlim(0,5)
        plt.plot(self.T,y)

        plt.ylabel('Specific heat per spin m')
        plt.xlabel('Temperature T')
        plt.show()

L55 = Ising(5)
L55.generator_s()
L55.calc()

y = L55.mag/max(L55.mag)

plt.style.use('seaborn-whitegrid')
plt.ylim(0,1.1)
plt.xlim(0,5)
plt.plot(L55.T,y)

plt.ylabel('Abosolute magnetization per spin m')
plt.xlabel('Temperature T')
plt.savefig("L55_mag.png")

plt.close()
y = L55.spec/max(L55.spec)

plt.style.use('seaborn-whitegrid')
plt.ylim(0,1.1)
plt.xlim(0,5)
plt.plot(L55.T,y)

plt.ylabel('Specific heat per spin m')
plt.xlabel('Temperature T')
plt.savefig("L55_spec.png")