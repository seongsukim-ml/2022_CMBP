import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.integrate as integrate
import mpmath as mp
import Exact_Ising_model
import os
import csv

Exact_M_data, Exact_C_data = Exact_Ising_model.Exact_calc()

def draw_fig1(path):
    a = pd.read_csv(path)
    T2 = a.iloc[:,1].values
    M2 = a.iloc[:,2].values
    C2 = a.iloc[:,3].values
    M2error = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))

    plt.style.use('seaborn-whitegrid')
    plt.figure(dpi=300)


    plt.ylim(-0.1,2)
    plt.xlim(0,T2[-1])
    # plt.errorbar(T2,M2,yerr=M2error, lw=1, linestyle='', marker='s', markersize=5,capsize=7, color='b', label='magnetization',mfc='none')
    plt.plot(T2,M2, lw=1, linestyle='', marker='s', markersize=5, color='b', label='magnetization',mfc='none')
    plt.plot(T2,C2,linestyle='', marker='o', markersize=5, color='orange', label='specific heat',mfc='none')
    # plt.errorbar(T2,C2,yerr=C2error,linestyle='', marker='o', markersize=5,capsize=7, color='orange' ,label='specific heat',mfc='none')


    plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    plt.legend()
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')


    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Temperature T')
    plt.show()

def draw_multi(path_list,Lsize,Bins,Step_Size):
    plt.style.use('seaborn-whitegrid')
    # plt.figure(dpi=300)

    for i, data in enumerate(path_list):
        a = pd.read_csv(data)
        T2 = a.iloc[:,1].values
        M2 = a.iloc[:,2].values
        C2 = a.iloc[:,3].values
        M2error = np.sqrt(1/(Step_Size[i]-1)*abs(a.iloc[:,5].values-(a.iloc[:,4].values)**2))
        # plt.ylim(-0.1,2)
        # plt.xlim(0,T2[-1])
        # plt.errorbar(T2,M2,yerr=M2error, lw=1, linestyle='', marker='s', markersize=5,capsize=7, color='b', label='magnetization',mfc='none')
        plt.errorbar(T2,M2,M2error,lw=1, linestyle='', marker='s', markersize=5, label='magnetization'+str(Lsize[i]),mfc='none')
        plt.plot(T2,C2,linestyle='', marker='o', markersize=5, label='specific heat'+str(Lsize[i]),mfc='none')
        # plt.errorbar(T2,C2,yerr=C2error,linestyle='', marker='o', markersize=5,capsize=7, color='orange' ,label='specific heat',mfc='none')

    plt.legend()
    plt.ylim(-0.1,3)
    # plt.xlim(T2[0],T2[-1])

    plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')


    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Temperature T')
    plt.show()

def draw_multi_with_exact(path_list,Lsize):
    plt.style.use('seaborn-whitegrid')
    # plt.figure(dpi=300)

    for i, data in enumerate(path_list):
        a = pd.read_csv(data)
        T2 = a.iloc[:,1].values
        M2 = a.iloc[:,2].values
        C2 = a.iloc[:,3].values
        M2error = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))

        # plt.ylim(-0.1,2)
        # plt.xlim(0,T2[-1])
        # plt.errorbar(T2,M2,yerr=M2error, lw=1, linestyle='', marker='s', markersize=5,capsize=7, color='b', label='magnetization',mfc='none')
        plt.plot(T2,M2, lw=1, linestyle='', marker='s', markersize=5, label='magnetization'+str(Lsize[i]),mfc='none')
        plt.plot(T2,C2,linestyle='', marker='o', markersize=5, label='specific heat'+str(Lsize[i]),mfc='none')
        # plt.errorbar(T2,C2,yerr=C2error,linestyle='', marker='o', markersize=5,capsize=7, color='orange' ,label='specific heat',mfc='none')

    plt.legend()
    plt.plot(*Exact_M_data, color='b')
    plt.plot(*Exact_C_data, color='orange')
    plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')


    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Temperature T')
    plt.show()

def draw_binder_FFS(path_list, Lsize, Bins, Step_Size, xliml = None, yliml = (0.75,1), Tc=2.269, nu=1):
    Llist = []
    for path in path_list:
        Llist.append(pd.read_csv(path))
    marker = ["o","s","v","^","8"]
    T2 = Llist[0].iloc[:,1].values
    Binder = [0.5*(3-i.iloc[:,6].values/(i.iloc[:,5].values)**2) for i in Llist]
    if xliml is None:
        xliml = (T2[0],T2[-1])

    plt.style.use('seaborn-whitegrid')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylim(yliml)
    plt.xlim(xliml)
    for i in range(len(path_list)):
        plt.plot(T2,Binder[i],marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i]))
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    plt.legend()
    plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    plt.ylabel('Binder ratio g')
    plt.xlabel('Temperature T')
    # plt.text(2.265,0.748, '$T_c$')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylim(0.75,1)
    plt.xlim(-2,2)
    T22 = T2-Tc
    for i in range(len(path_list)):
        plt.plot(T22*(Lsize[i]**(1/nu)),Binder[i],linestyle="",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i]))
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.legend()
    plt.axvline(x=0,c='grey',lw=1,dashes=[2,2])
    plt.ylabel('Binder ratio g')
    plt.xlabel('Finite size scailing')
    # plt.xlabel('$L^{1/\nu}[T-T_c]$')
    plt.show()

def improved_list_dir(path, rms = ["Input"]):
    ls = os.listdir(path)
    ls.sort()

    for rm in rms:
        if rm in ls:
            ls.pop(ls.index(rm))
    a100 = "a=100"
    if a100 in ls:
        ls.append(ls.pop(ls.index(a100)))
    return ls

def parsing_helper(path):
    parsing_dict = {}
    alpha_list = []
    for dir in improved_list_dir(path):
        ls = os.listdir(path+dir+'/')
        Lsize, Bins, mcs = [],[],[]
        ls.sort()

        for i, file in enumerate(ls):
            comp = file.split('_')
            Lsize.append(int(comp[3]))
            Bins.append(int(comp[5][3:]))
            mcs.append(int(comp[6][3:]))
            ls[i] = path+dir+'/' + file
        alpha_list.append(dir)
        parsing_dict[dir]= [ls, Lsize, Bins, mcs]
    return parsing_dict, alpha_list

crit_range_dict = {"a=1"    :[1.58,1.64],
                   "a=1.5"  :[1.76,1.82],
                   "a=2"    :[1.88,1.94],
                   "a=2.5"  :[1.96,2.02],
                   "a=3"    :[2.06,2.12],
                   "a=3.5"  :[2.08,2.14],
                   "a=100"  :[2.22,2.28],
}

def draw_1binder(path):
    L5 = pd.read_csv(path)
    # L5 = pd.read_csv("../C++/Result/Exact_c_5_int40_1.csv")
    Llist = [L5]
    T2 = L5.iloc[:,1].values
    Binder = [0.5*(3-i.iloc[:,6].values/(i.iloc[:,5].values)**2) for i in Llist]
    plt.style.use('seaborn-whitegrid')

    # plt.ylim(-0.1,2)
    plt.xlim(1.5,4)
    marker = ["o","s","o","s","o"]
    for i in range(1):
        plt.plot(T2,Binder[i],marker=marker[i],markersize=5)

    # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])

    plt.show()