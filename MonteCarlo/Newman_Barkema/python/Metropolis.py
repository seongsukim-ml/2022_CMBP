import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd

## Need to be fixed

class Metropolis:
    ## In 2D Ising model => 2/log(1+sqrt(2))
    T_CRIT = 2.269185

    def __init__(self, L, Bin, Tsrt, Tfin, B = 0, J =1, isTinf=False):
        self.L = L
        self.N = L*L
        self.Bin = Bin
        self.B = B
        self.J = J
        self.isTinf = isTinf

        self.XNN = 1
        self.YNN = L

        self.HH = 0
        self.sigma = 0

        self.MV = np.array(Bin,np.double)
        self.CV = np.array(Bin,np.double)
        self.TV = np.array(Bin,np.double)
        self.BetaV = np.array(Bin,np.double)

        self.sc = np.ones(self.N,dtype=np.int0)
        self.prob = np.zeros(5,dtype=np.double)

        self.Fliped_Step = 0
        self.Total_Step = 0
        self.Calc_call = 0

        for i in range(Bin):
            if Bin == 1:
                self.TV = Tsrt
                self.BetaV = 1/self.TV
                break
            elif Tsrt != 0:
                self.TV[i] = Tsrt + ((Tfin-Tsrt)/(Bin))*(i+1)
            else :
                self.TV[i] = Tsrt + ((Tfin-Tsrt)/(Bin-1))*i
            self.BetaV[i] = 1/self.TV[i]
    
    def Prob_calc(self,beta):
        for i in range(2):
            self.prob[2*(i+1)] = np.exp(-2*beta*2*(i+1))

    def Initialize(self,beta):
        self.Calc_call = 0
        self.Fliped_Step = 0
        self.Total_Step = 0

        self.sc = np.ones(self.N,dtype=np.int0)
        if(self.isTinf):
            self.sc = np.array([1-(int)(np.random.rand()*2)*2 for i in range(self.N)],dtype=np.int0)

        self.Prob_calc()
        self.Measure()

    def SweepHelical(self,i):
        sum = 0

        nn = i - 1
        if(nn < 0) : nn += self.N
        sum += self.sc[nn]

        nn = i + 1
        if(nn >= self.N): nn -= self.N
        sum += self.sc[nn]

        nn = i - self.L
        if(nn < 0): nn += self.N
        sum += self.sc[nn]

        nn = i + self.L
        if(nn >= self.N): nn -= self.N
        sum += self.sc[nn]
        return sum

    def BoundaryHelical(self,i):
        sum = 0

        nn = i + 1
        if(nn == self.N): nn = 0
        sum += self.sc[nn]

        nn = i + self.L
        if(nn >= self.N): nn -= self.N
        sum += self.sc[nn]
        return sum

    def SweepPBC(self,i):
        sum = 0

        nn = i -1
        if((nn+1 % self.L) == 0) : nn += self.L
        sum += self.sc[nn]

        nn = i + 1
        if(nn % self.L == 0): nn -= self.L
        sum += self.sc[nn]

        nn = i - self.L
        if(nn < 0): nn += self.N
        sum += self.sc[nn]

        nn = i + self.L
        if(nn >= self.N): nn -= self.N
        sum += self.sc[nn]
        return sum

    def Measure(self,func = BoundaryHelical):
        res = 0
        for i in range(self.N):
            sum = func(self,i)
            res += self.J*sum*self.sc[i]

        sigma = np.sum(self.sc)

        HH = -res -self.B*sigma
        self.HH = HH
        self.sigma = sigma

        return HH, sigma

    def MeasureFast(self):
        return self.HH, self.sigma    

    def Calculate(self,func = SweepHelical, Random=False):
        for i in range(self.N):
            if(Random):
                k = np.random.randint(0,self.N-1)
            elif (self.N%2 ==0):
                k = 2*i
                if(k < self.N): k = k if int(k/self.L)%2 ==0 else k+1
                else:
                    k -= self.N
                    k = k+1 if int(k/self.L)%2 ==0 else k
            else:
                k = 2*i if 2*i < self.N else 2*i-self.N

            # delta = Enew - Eold
            # delta = 2*self.sc[k]*self.J*self.SweepHelical(k)
            # delta = self.J*func(k)*self.sc[k]
            delta = func(self,k)*self.sc[k]
            self.Total_Step += 1

            # print(delta)
            if(delta <= 0 or (np.random.rand() < self.prob[delta])): # A = 1
                self.Fliped_Step += 1
                self.sc[k] *= -1
                self.sigma += 2*self.sc[k]
                self.HH += 2*delta

    def IteraterUntilEquilibrium(self, equil_time, func = SweepHelical, Random =False):
        for i in range(equil_time):
            self.Calculate(func = func, Random=Random)