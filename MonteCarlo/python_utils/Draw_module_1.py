from tkinter import BOTTOM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.integrate as integrate
import mpmath as mp
import Exact_Ising_model
import os
import csv
import copy

from scipy.optimize import curve_fit

Exact_M_data, Exact_C_data = Exact_Ising_model.Exact_calc()

alpha_file_name_dict = {
                    "a=1"    :"1.000000",
                    "a=1.5"  :"1.500000",
                    "a=2"    :"2.000000",
                    "a=2.5"  :"2.500000",
                    "a=3"    :"3.000000",
                    "a=3.5"  :"3.500000",
                    "a=100"  :"100.000000",
}

plt.style.use('seaborn-whitegrid')
hr = {"figsize":(10,6),"dpi":500}

def draw_fig1(path, Lsize, plot_list = ["MM","CC"]):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if(isinstance(path,str)):
        a = pd.read_csv(path)
    else:
        a = path
    TT = a.iloc[:,1].values
    MM = a.loc[:,'abs(mm)'].values
    CC = a.loc[:,'specific heat'].values
    Sus = ((a.loc[:,'mm**2'].values)*Lsize**4 - (MM*Lsize**2)**2)
    MMerror = a.loc[:,'MMerr'].values
    CCerror = a.loc[:,'CCerr'].values
    # MMerror = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))
    Suserr  = abs(a.loc[:,'MM2err'].values**(0.5)*(Lsize**4)+abs(MM*2*MMerror*Lsize**4))
    # print(MM*2*MMerror)

    plt.style.use('seaborn-whitegrid')
    # plt.figure(dpi=300)


    # plt.ylim(-0.1,2)
    # plt.xlim(0,TT[-1])
    if "MM"    in plot_list:
        plt.plot(TT,MM, lw=1, linestyle='', marker='s', markersize=5, color='b', label='magnetization',mfc='none')
    if "CC"    in plot_list:
        plt.plot(TT,CC,linestyle='', marker='o', markersize=5, color='orange', label='specific heat',mfc='none')
    if "Sus"   in plot_list:
        plt.plot(TT,Sus,linestyle='', marker='o', markersize=5, color='r', label='mag susceptibility',mfc='none')
    if "MMerr" in plot_list:
        plt.errorbar(TT,MM,yerr=MMerror, lw=1, linestyle='', marker='s', markersize=5,capsize=7, color='b', label='magnetization',mfc='none')
    if "CCerr" in plot_list:
        plt.errorbar(TT,CC,yerr=CCerror,linestyle='', marker='o', markersize=5,capsize=7, color='orange' ,label='specific heat',mfc='none')
    if "Suserr"   in plot_list:
        plt.errorbar(TT,Sus,yerr=Suserr,linestyle='', marker='o', markersize=5,capsize=7, color='r', label='mag susceptibility',mfc='none')

    # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.legend()
    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Temperature T')
    plt.show()


def draw_multi(path_list,Lsize,Bins,Step_Size, plot_list = ["MM","CC"]):
    plt.style.use('seaborn-whitegrid')
    # plt.figure(dpi=300)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i, data in enumerate(path_list):
        if(isinstance(data,str)):
            a = pd.read_csv(data)
        else:
            a = data

        TT = a.iloc[:,1].values
        MM = a.loc[:,'abs(mm)'].values
        CC = a.loc[:,'specific heat'].values
        # Sus = (a.loc[:,'mm**2'].values) - MM**2
        Sus = (TT*Lsize[i]**2*(a.loc[:,'mm**2'].values) - (MM)**2)
        MMerror = a.loc[:,'MMerr'].values
        CCerror = a.loc[:,'CCerr'].values
        # MMerror = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))
        Suserr  = abs(a.loc[:,'MM2err'].values**(0.5)*(Lsize[i]**4)+abs(MM*2*MMerror*Lsize[i]**4))
        # print(MM*2*MMerror)

        plt.style.use('seaborn-whitegrid')
        # plt.figure(dpi=300)


        # plt.ylim(-0.1,2)
        # plt.xlim(0,TT[-1])
        if "MM"    in plot_list:
            plt.plot(TT,MM, lw=1, linestyle='', marker='s', markersize=5, label='magnetization'+str(Lsize[i]),mfc='none')
        if "CC"    in plot_list:
            plt.plot(TT,CC,linestyle='', marker='o', markersize=5, label='specific heat'+str(Lsize[i]),mfc='none')
        if "Sus"   in plot_list:
            plt.plot(TT,Sus,linestyle='', marker='o', markersize=5, label='mag susceptibility'+str(Lsize[i]),mfc='none')
        if "MMerr" in plot_list:
            plt.errorbar(TT,MM,yerr=MMerror, lw=1, linestyle='', marker='s', markersize=5,capsize=7, label='magnetization'+str(Lsize[i]),mfc='none')
        if "CCerr" in plot_list:
            plt.errorbar(TT,CC,yerr=CCerror,linestyle='', marker='o', markersize=5,capsize=7 ,label='specific heat'+str(Lsize[i]),mfc='none')
        if "Suserr"   in plot_list:
            plt.errorbar(TT,Sus,yerr=Suserr,linestyle='', marker='o', markersize=5,capsize=7, label='mag susceptibility'+str(Lsize[i]),mfc='none')

    if "exact" in plot_list:
        plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Temperature T')
    plt.legend()
    plt.show()

def draw_multi_critical(path_list,Lsize,Bins,Step_Size,Tc = 2.269, plot_list = ["MM","CC"]):
    plt.style.use('seaborn-whitegrid')
    # plt.figure(dpi=300)
    fig = plt.figure(figsize=(10,6),dpi=500)
    ax = fig.add_subplot(111)

    for i, data in enumerate(path_list):
        if(isinstance(data,str)):
            a = pd.read_csv(data)
        else:
            a = data

        TT = a.iloc[:,1].values
        MM = a.loc[:,'abs(mm)'].values
        # CC = a.loc[:,'specific heat'].values/Lsize[i]**(0.5)
        CC = a.loc[:,'specific heat'].values/np.log(Lsize[i]**2)

        # Sus = ((a.loc[:,'mm**2'].values)*Lsize[i]**4 - (MM*Lsize[i]**2)**2)
        # Sus = ((a.loc[:,'mm**2'].values)*Lsize[i]**4 - (MM*Lsize[i]**2)**2)*Lsize[i]**(-1.867*2)
        # Sus = ((a.loc[:,'mm**2'].values)*Lsize[i]**4 - (MM*Lsize[i]**2)**2)*Lsize[i]**(-1.75*2)
        # Sus = (TT*Lsize[i]**2*(a.loc[:,'mm**2'].values) - (MM)**2)*Lsize[i]**(-1.75)
        Sus = (TT*Lsize[i]**2*(a.loc[:,'mm**2'].values) - (MM)**2)*Lsize[i]**(-1.75)
        


        MMerror = a.loc[:,'MMerr'].values
        # CCerror = a.loc[:,'CCerr'].values/Lsize[i]**(0.5)
        CCerror = a.loc[:,'CCerr'].values/np.log(Lsize[i]**2)

        # MMerror = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))
        Suserr  = abs(a.loc[:,'MM2err'].values**(0.5)*(Lsize[i]**4)+abs(MM*2*MMerror*Lsize[i]**4))
        # print(MM*2*MMerror)

        plt.style.use('seaborn-whitegrid')
        # plt.figure(dpi=300)

        TT = (TT-Tc)*Lsize[i]/TT
        # plt.ylim(-0.1,2)
        # plt.xlim(0,TT[-1])
        if "MM"    in plot_list:
            plt.plot(TT,MM, lw=1, linestyle='', marker='s', markersize=5, label='magnetization'+str(Lsize[i]),mfc='none')
        if "CC"    in plot_list:
            plt.plot(TT,CC,linestyle='', marker='o', markersize=5, label='specific heat'+str(Lsize[i]),mfc='none')
        if "Sus"   in plot_list:
            plt.plot(TT,Sus,linestyle='', marker='o', markersize=5, label='mag susceptibility'+str(Lsize[i]),mfc='none')
        if "MMerr" in plot_list:
            plt.errorbar(TT,MM,yerr=MMerror, lw=1, linestyle='', marker='s', markersize=5,capsize=7, label='magnetization'+str(Lsize[i]),mfc='none')
        if "CCerr" in plot_list:
            plt.errorbar(TT,CC,yerr=CCerror,linestyle='', marker='o', markersize=5,capsize=7 ,label='specific heat'+str(Lsize[i]),mfc='none')
        if "Suserr"   in plot_list:
            plt.errorbar(TT,Sus,yerr=Suserr,linestyle='', marker='o', markersize=5,capsize=7, label='mag susceptibility'+str(Lsize[i]),mfc='none')

    if "exact" in plot_list:
        plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    # plt.text(0.05,1.91,"100 $\\times$ 100 lattice")
    # plt.text(0.05,1.79,"15000 steps per site")
    # plt.text(0.05,1.67,"Jackknife bin : 50")
    # plt.text(2.27, -0.2, '$T_c$', ha='center')

    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('Specific heat or magnetization per spin m')
    plt.xlabel('Reduced beta')
    plt.legend()
    plt.show()

def draw_binder_FFS(path_list, Lsize, Bins, Step_Size, xliml = None, yliml = (0.75,1), Tc=2.269, nu=1, _title = ""):
    Llist = []
    for data in path_list:
        if(isinstance(data,str)):
            Llist.append(pd.read_csv(data))
        else:
            Llist.append(data)
    marker = ["o","s","v","^","8"]
    TT = Llist[0].iloc[:,1].values
    Binder = [0.5*(3-i.iloc[:,6].values/(i.iloc[:,5].values)**2) for i in Llist]
    if xliml is None:
        xliml = (TT[0],TT[-1])

    plt.style.use('seaborn-whitegrid')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylim(yliml)
    plt.xlim(xliml)
    for i in range(len(path_list)):
        plt.plot(TT,Binder[i],marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i]))
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
    T22 = TT-Tc
    for i in range(len(path_list)):
        plt.plot(T22*(Lsize[i]**(1/nu)),Binder[i],linestyle="",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i]))
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.legend()
    plt.axvline(x=0,c='grey',lw=1,dashes=[2,2])
    plt.ylabel('Binder ratio g')
    plt.xlabel('Finite size scailing')
    plt.title(_title)
    # plt.xlabel('$L^{1/\nu}[T-T_c]$')
    plt.show()

def draw_binder_multi(path_list,Lsize,Bins,Step_Size, Tc=2.269, nu=1 , plot_list = ["MM","CC","hr"],xliml=[],yliml = []):
    if "hr" in plot_list:
        fig = plt.figure(**hr)
    else:
        fig = plt.figure()

    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o"]

    for i,path in enumerate(path_list):
        if(isinstance(path,str)):
            a = pd.read_csv(path)
        else:
            a = path

        TT = a.iloc[:,1].values
        if "calc" in plot_list:
            Binder = 0.5*(3-a.loc[:,"mm**4"].values/(a.loc[:,"mm**2"].values)**2)*(2/3)
        else:
            Binder = a.loc[:,'Binder'].values*(2/3)
        # plt.ylim(-0.1,2)
        # plt.xlim(1.5,4)
        plt.style.use('seaborn-whitegrid')
        if "BBerr" in plot_list:
            BBerr = a.loc[:,"BBerr"].values*(2/3)
            plt.errorbar(TT,Binder,BBerr,marker=marker[0],mfc='none',markersize=5,capsize= 7,label = "L"+str(Lsize[i]))
        else :
            plt.plot(TT,Binder,marker=marker[0],mfc='none',markersize=5)
        # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    if xliml:
        plt.xlim(xliml)
    if yliml:
        plt.ylim(yliml)
    plt.legend()
    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('Binder ratio g')
    plt.xlabel('Temperature')

    plt.show()

def draw_binder_FFS_multi(path_list,Lsize,Bins,Step_Size, Tc=2.269, nu=1 , plot_list = ["MM","CC"]):
    fig = plt.figure(figsize=(10,6),dpi=500)
    # fig = plt.figure()
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o"]

    for i,path in enumerate(path_list):
        if(isinstance(path,str)):
            a = pd.read_csv(path)
        else:
            a = path

        TT = a.iloc[:,1].values
        if "calc" in plot_list:
            print("calc")
            Binder = 0.5*(3-a.loc[:,"mm**4"].values/(a.loc[:,"mm**2"].values)**2)*(2/3)
        else:
            Binder = a.loc[:,'Binder'].values*(2/3)
        # plt.ylim(-0.1,2)
        # plt.xlim(1.5,4)
        plt.style.use('seaborn-whitegrid')
        TTa = TT-Tc

        if "BBerr" in plot_list:
            BBerr = a.loc[:,"BBerr"].values*(2/3)
            plt.errorbar(TTa*(Lsize[i]**(1/nu)),Binder,BBerr,linestyle="",marker=marker[0],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i]))
        else :
            plt.plot(TTa*(Lsize[i]**(1/nu)),Binder,linestyle="",marker=marker[0],markersize=5,mfc='none',label="L"+str(Lsize[i]))
        # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])

    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('Binder ratio g')
    plt.xlabel('(T-Tc)*(Lsize^(1/nu))')
    plt.show()

def draw_binder_fig1(path, plot_list = ["BBerr"]):
    # options : "BBerr", "calc"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if(isinstance(path,str)):
        a = pd.read_csv(path)
    else:
        a = path

    TT = a.iloc[:,1].values
    if "calc" in plot_list:
        Binder = 0.5*(3-a.loc[:,"mm**4"].values/(a.loc[:,"mm**2"].values)**2)
    else:
        Binder = a.loc[:,'Binder'].values
    # plt.ylim(-0.1,2)
    # plt.xlim(1.5,4)
    marker = ["o","s","o","s","o"]
    plt.style.use('seaborn-whitegrid')
    if "BBerr" in plot_list:
        BBerr = a.loc[:,"BBerr"].values
        plt.errorbar(TT,Binder,BBerr,marker=marker[0],mfc='none',markersize=5,capsize= 7)
    else :
        plt.plot(TT,Binder,marker=marker[0],mfc='none',markersize=5)
    # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.show()

def draw_binder_FSS_fig1(path,Tc, Lsize,nu = 1,plot_list = ["BBerr"]):
    # options : "BBerr", "calc"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if(isinstance(path,str)):
        a = pd.read_csv(path)
    else:
        a = path

    TT = a.iloc[:,1].values
    if "calc" in plot_list:
        Binder = 0.5*(3-a.loc[:,"mm**4"].values/(a.loc[:,"mm**2"].values)**2)
    else:
        Binder = a.loc[:,'Binder'].values
    # plt.ylim(-0.1,2)
    # plt.xlim(1.5,4)
    marker = ["o","s","o","s","o"]
    plt.style.use('seaborn-whitegrid')
    TTa = TT-Tc

    if "BBerr" in plot_list:
        BBerr = a.loc[:,"BBerr"].values
        # plt.errorbar(TT,Binder,BBerr,marker=marker[0],markersize=5, capsize= 7)
        plt.errorbar(TTa*(Lsize**(1/nu)),Binder,BBerr,linestyle="",marker=marker[0],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize))

    else :
        # plt.plot(TT,Binder,marker=marker[0],markersize=5)
        plt.plot(TTa*(Lsize**(1/nu)),Binder,linestyle="",marker=marker[0],markersize=5,mfc='none',label="L"+str(Lsize))

    # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.show()

def draw_correlation(path_list,Lsize,Bins,Step_Size,err_data = [],idx_except=[],plot_list = ['long','short','err','abs','except0','half']):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o","s","o","s"]*2

    for i,path in enumerate(path_list):
        if i in idx_except:
            continue
        cur_L = Lsize[i]
        cur_data = copy.deepcopy(path)
        if 'err' in plot_list:
            longerr  = err_data[i][0]
            shorterr = err_data[i][1]
        if 'half' in plot_list:
            cur_L = int(cur_L/2)
            cur_data[0] = cur_data[0][:cur_L]
            cur_data[1] = cur_data[1][:cur_L]
            if 'err' in plot_list:
                longerr  = longerr[:cur_L]
                shorterr = shorterr[:cur_L]
        Lrange = range(cur_L)
        if 'except0' in plot_list:
            cur_L -= 1
            cur_data[0] = cur_data[0][1:]
            cur_data[1] = cur_data[1][1:]
            if 'err' in plot_list:
                longerr  = longerr[1:]
                shorterr = shorterr[1:]
            Lrange = range(1,cur_L+1)
        if 'long' in plot_list:
            if 'abs' in plot_list:
                cur_data[0] = abs(cur_data[0])
            if 'err' in plot_list:
                print(len(cur_data[0]),cur_L,len(longerr))
                plt.errorbar(Lrange,cur_data[0],longerr,linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
            else:
                plt.plot(Lrange,cur_data[0],linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
        if 'short' in plot_list:
            if 'err' in plot_list:
                plt.errorbar(Lrange,cur_data[1],shorterr,linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
            else:
                plt.plot(Lrange,cur_data[1],linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')

    plt.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    # plt.ylabel('Binder ratio g')
    # plt.xlabel('(T-Tc)*(Lsize^(1/nu))')
    plt.show()

def draw_correlation_advanced_1(path_list,Lsize,Bins,Step_Size,err_data = [],idx_except=[],plot_list = ['long','short','err','abs','except0','half','log_fit'],temp = 0):
    fig = plt.figure(figsize=(10,6),dpi=500)
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o","s","o","s"]*2

    for i,path in enumerate(path_list):
        if i in idx_except:
            continue
        cur_L = Lsize[i]
        cur_data = copy.deepcopy(path)
        if 'abs' in plot_list:
            cur_data[0] = abs(cur_data[0])
        if 'err' in plot_list:
            longerr  = err_data[i][0]
            shorterr = err_data[i][1]
        if 'half' in plot_list:
            cur_L = int(cur_L/2)+1
            cur_data[0] = cur_data[0][:cur_L]
            cur_data[1] = cur_data[1][:cur_L]
            if 'err' in plot_list:
                longerr  = longerr[:cur_L]
                shorterr = shorterr[:cur_L]
        Lrange = range(cur_L)
        if 'except0' in plot_list:
            cur_L -= 1
            cur_data[0] = cur_data[0][1:]
            cur_data[1] = cur_data[1][1:]
            if 'err' in plot_list:
                longerr  = longerr[1:]
                shorterr = shorterr[1:]
            Lrange = range(1,cur_L+1)
        if 'log' in plot_list:
            if 'err' in plot_list:
                longerr  = abs(longerr/np.log(10)/cur_data)
                shorterr = abs(shorterr/np.log(10)/cur_data)
            Lrange = np.log10(Lrange)
            cur_data = np.log10(cur_data)
        if 'log_fit' in plot_list:
            ff = lambda X,a,b:X**(-a)+b
            print(np.array(Lrange),cur_data)
            popt, popv = curve_fit(ff,np.array(Lrange),cur_data[0])
        # print(cur_data[0],longerr)
        if 'long' in plot_list:
            if 'dim_scale' in plot_list:
                if 'log_s' in plot_list:
                    plt.errorbar(Lrange/Lsize[i],(cur_data[0]*Lsize[i]**(0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
                    # plt.plot(Lrange/Lsize[i],(cur_data[0]*Lsize[i]**(0.25)),linestyle=" ",marker=marker[i],mfc='none',label="L"+str(Lsize[i])+'l')
                    # plt.ylim((10**-3,1.5))
                    plt.yscale('log')
                    plt.ylabel('$C(r)L^{1/4}$')
                    plt.xlabel('r/L')
                else:
                    plt.errorbar(Lrange/Lsize[i],(cur_data[0]*Lsize[i]**(0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
                    plt.ylabel('$C(r)L^{1/4}$')
                    plt.xlabel('r/L')
            elif 'cft' in plot_list:
                plt.errorbar(Lrange/Lsize[i],cur_data[0]*(Lsize[i]*(np.sin(np.pi*np.array(Lrange/Lsize[i]))))**(0.25),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=4,mfc='none',label="L"+str(Lsize[i])+'l')
                # plt.errorbar(Lrange/Lsize[i],(cur_data[0])/(np.sin(3.141592*np.array(Lrange/Lsize[i]))**(-0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=4,mfc='none',label="L"+str(Lsize[i])+'l')
                plt.ylabel('$C(r)[L\sin(\pi r/L)]^{1/4}$')
                plt.xlabel('r/L')
            elif 'exponent' in plot_list:
                plt.errorbar(Lrange,cur_data[0]/(np.array(Lrange))**(-0.25),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=4,mfc='none',label="L"+str(Lsize[i])+'l')
            elif 'err' in plot_list:
                print(len(cur_data[0]),cur_L,len(longerr))
                plt.plot(Lrange/Lsize[i],cur_data[0],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'l')
            elif 'loglog' in plot_list:
                plt.loglog(Lrange,cur_data[0],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'l')
            else:
                plt.plot(Lrange,cur_data[0],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'l')
            if 'log_fit' in plot_list:
                plt.loglog(Lrange,ff(Lrange,*popt))
                print(popt)
            # if 'cft_test' in plot_list:
                # plt.errorbar(np.sin(Lrange/Lsize[i]),cur_data[0]*Lsize[i]**(0.25),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')

        if 'short' in plot_list:
            if 'dim_scale' in plot_list:
                plt.errorbar(Lrange/Lsize[i],(cur_data[1]*Lsize[i]**(0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
            elif 'cft' in plot_list:
                plt.errorbar(Lrange/Lsize[i],(cur_data[1]*Lsize[i]**(0.25))/(np.sin(3.141592*np.array(Lrange/Lsize[i]))**(-0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
            elif 'err' in plot_list:
                print(len(cur_data[0]),cur_L,len(longerr))
                plt.errorbar(Lrange/Lsize[i],(cur_data[1]*Lsize[i]**(0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
                # plt.errorbar(Lrange/Lsize[i],(cur_data[1]*Lsize[i]**(0.25))/(np.sin(3.141592*np.array(Lrange/Lsize[i]))**(-0.25)),longerr,linestyle=" ",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
                # plt.plot(np.sin(3.141592*np.array(Lrange/Lsize[i]))**(-0.25),cur_data[1]*Lsize[i]**(0.25),linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'l')
            # if 'err' in plot_list:
                # plt.errorbar(Lrange,cur_data[1],shorterr,linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
            else:
                plt.plot(Lrange,cur_data[1],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'s')

    plt.legend()
    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    # plt.ylabel('Binder ratio g')
    # plt.xlabel('(T-Tc)*(Lsize^(1/nu))')
    plt.show()

def draw_correlation_advanced(path_list,Lsize,Bins,Step_Size,err_data = [],idx_except=[],plot_list = ['long','short','err','abs','except0','half']):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o","s","o","s"]*2

    for i,path in enumerate(path_list):
        if i in idx_except:
            continue
        cur_L = Lsize[i]
        cur_data = copy.deepcopy(path)
        if 'abs' in plot_list:
            cur_data[0] = abs(cur_data[0])
        if 'err' in plot_list:
            longerr  = err_data[i][0]
            shorterr = err_data[i][1]
        if 'half' in plot_list:
            cur_L = int(cur_L/2)
            cur_data[0] = cur_data[0][:cur_L]
            cur_data[1] = cur_data[1][:cur_L]
            if 'err' in plot_list:
                longerr  = longerr[:cur_L]
                shorterr = shorterr[:cur_L]
        Lrange = range(cur_L)
        temp = range(cur_L)
        if 'except0' in plot_list:
            cur_L -= 1
            cur_data[0] = cur_data[0][1:]
            cur_data[1] = cur_data[1][1:]
            if 'err' in plot_list:
                longerr  = longerr[1:]
                shorterr = shorterr[1:]
            Lrange = range(1,cur_L+1)
            temp = range(1,cur_L+1)
        Lrange = np.log(Lrange)
        # print(Lrange)
        # print(cur_data)
        cur_data = np.log(cur_data)
        temp1 = np.array(range(4,cur_L,4))
        temp2 = np.array(range(1,len(temp1)+1))
        # print(temp1,temp2)
        # print(len(temp1),len(temp2))
        # print(cur_data[0][temp2])
        new_data = [[],[]]
        new_data[0] = cur_data[0][temp2]
        new_data[1] = cur_data[1][temp2]
        print(new_data)
        cur_data = new_data
        Lrange = np.log(temp1)
        if 'err' in plot_list:
            longerr  = np.log(longerr)
            shorterr = np.log(shorterr)
        if 'long' in plot_list:
            if 'err' in plot_list:
                print(len(cur_data[0]),cur_L,len(longerr))
                plt.errorbar(Lrange,cur_data[0],longerr,linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'l')
            else:
                plt.plot(Lrange,cur_data[0],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'l')
        if 'short' in plot_list:
            if 'err' in plot_list:
                plt.errorbar(Lrange,cur_data[1],shorterr,linestyle="-",marker=marker[i],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i])+'s')
            else:
                plt.plot(Lrange,cur_data[1],linestyle="-",marker=marker[i],markersize=5,mfc='none',label="L"+str(Lsize[i])+'s')

    plt.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('$ln|C(L/4)|$')
    plt.xlabel('$ln L$')
    plt.show()

def draw_correlation_fitting(path_list,Lsize,Bins,Step_Size,err_data = [],idx_except=[],plot_list = ['long','short','err','abs','except0','half']):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o","s","o","s"]*2
    xx,yy,err = [], [], []
    for i,path in enumerate(path_list):
        if i in idx_except:
            continue
        cur_L = int(Lsize[i]/8)
        cur_data = copy.deepcopy(path)
        longerr  = err_data[i][0][cur_L]
        xx.append(np.log(cur_L))
        yy.append(np.log(abs(cur_data[0][cur_L])))
        err.append(longerr/cur_data[0][cur_L])
    xx = np.array(xx)
    yy = np.array(yy)
    err= np.array(err)
    ff = lambda x,a,b: a*x+b
    popt,popv = curve_fit(ff,xx,yy)
    print(err)
    plt.errorbar(xx,yy,yerr=err,label='data',linestyle=" ",marker='o',markersize=5,capsize=7)
    plt.plot(xx,ff(xx,*popt),label = 'fitting')
    print(popt)
    print(np.sqrt(popv))

    plt.legend()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('$ln|C(L/4)|$')
    plt.xlabel('$ln L$')
    plt.show()


######################################

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
        ls = sorted(ls,key=cmp_to_key(file_name_compare))
        # print(ls)

        for i, file in enumerate(ls):
            comp = file.split('_')
            Lsize.append(int(comp[3]))
            Bins.append(int(comp[5][3:]))
            mcs.append(int(comp[6][3:]))
            ls[i] = pd.read_csv(path+dir+'/' + file)
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

from functools import cmp_to_key

def file_name_compare(a,b):
    asp = a.split('_')
    bsp = b.split('_')


    stdidx = 3
    cnt = 0
    for aa in asp:
        if "int" in aa:
            stdidx = cnt
        cnt += 1
    # print(stdidx)
    # stdidx = asp.find('int') ## binning element, ex) int21
    if(int(asp[stdidx-1]) > int(bsp[stdidx-1])):
        # print(asp[stdidx-1],bsp[stdidx-1])
        return 1
    elif int(asp[stdidx-1]) == int(bsp[stdidx-1]):
        if asp[stdidx] > bsp[stdidx] :
            return 1
        elif asp[stdidx] == bsp[stdidx]:
            if asp[stdidx+1] > bsp[stdidx+1]:
                return 1
            elif asp[stdidx+1] == bsp[stdidx+1]:
                asp[-1] = int(asp[-1][:-4])
                bsp[-1] = int(bsp[-1][:-4])
                if asp[-1] > bsp[-1]:
                    return 1
                elif asp[-1] == bsp[-1]:
                    return 0
                else:
                    return -1
            else:
                return -1
        else:
            return -1
    else:
        return -1


def parsing_helper_cor(path):
    parsing_dict = {}
    parsing_dict_corres = {}
    parsing_dict_corerr = {}
    alpha_list = []
    for dir in improved_list_dir(path):
        ls = os.listdir(path+dir+'/')
        Lsize, Bins, mcs = [],[],[]
        ls = sorted(ls,key=cmp_to_key(file_name_compare))
        # print(ls)

        ls_corerr, Lsize_corerr, Bins_corerr, mcs_corerr = [], [], [], []
        ls_corres, Lsize_corres, Bins_corres, mcs_corres = [], [], [], []
        ls_std    = []
        for i, file in enumerate(ls):
            comp = file.split('_')
            res = pd.read_csv(path+dir+'/' + file)
            if(file.find('corres') != -1):
                Lsize_corres.append(int(comp[3]))
                Bins_corres.append(int(comp[5][3:]))
                mcs_corres.append(int(comp[6][3:]))
                ls_corres.append(res)
            elif(file.find('corerr') != -1):
                Lsize_corerr.append(int(comp[3]))
                Bins_corerr.append(int(comp[5][3:]))
                mcs_corerr.append(int(comp[6][3:]))
                ls_corerr.append(res)
            else:
                Lsize.append(int(comp[3]))
                Bins.append(int(comp[5][3:]))
                mcs.append(int(comp[6][3:]))
                ls_std.append(res)

        alpha_list.append(dir)
        # print(ls)
        parsing_dict_corerr[dir] = [ls_corerr, Lsize_corerr, Bins_corerr, mcs_corerr]
        parsing_dict_corres[dir] = [ls_corres, Lsize_corres, Bins_corres, mcs_corres]
        parsing_dict[dir] = [ls_std,Lsize,Bins,mcs]
    return parsing_dict, parsing_dict_corres, parsing_dict_corerr, alpha_list

def cor_data_generator(parsing_dict_cor,alpha_list,option = []):
    new_parsing_dict_cor = {}
    for alpha in alpha_list:
        a = parsing_dict_cor[alpha]
        ls, Lsize, Bins, mcs = a[0], a[1], a[2], a[3]
        Luniq = pd.array(Lsize).unique()
        idx = 0
        res = [[],[],[],[]]
        # print(idx)
        for ll in Luniq:
            cnt = Lsize.count(ll)
            res[1].append(ll)
            long_ls,short_ls = [], []
            for ii in range(idx,idx+cnt):
                # print(ls[ii])
                long, short = cor_build(ls[ii],ll,option)
                long_ls.append(long)
                short_ls.append(short)
            tt = [merge_upti_data(long_ls,option),merge_upti_data(short_ls,option)]
            res[0].append(tt)
            res[2].append(Bins[0])
            res[3].append(sum(mcs[idx:idx+cnt]))
            idx += cnt
        new_parsing_dict_cor[alpha] = res
    return new_parsing_dict_cor


def cor_build(data,Lsize,option = []):
    long_res = []
    short_res = []
    for rr in range(data.shape[0]):
        if 'full' in option:
            restore = restore_1Dvec_to_matrix
        else:
            restore = restore_upper_triangle
        xLsize = Lsize
        if '2Ly' in option:
            yLsize = 2*Lsize
        if '2Lx' in option:
            yLsize = Lsize//2
        if(rr%2 == 0):
            upti = restore(data.iloc[rr][2:],yLsize)
            long_res.append(upper_triangle_statics(upti,yLsize,option))
        else:
            upti = restore(data.iloc[rr][2:],xLsize)
            short_res.append(upper_triangle_statics(upti,xLsize,option))
    # return 1, 2
    return merge_upti_data(long_res,option), merge_upti_data(short_res,option)


# def cor_build(data,Lsize,option = []):
#     long_res = []
#     short_res = []
#     for rr in range(data.shape[0]):
#         if 'full' in option:
#             upti = restore_1Dvec_to_matrix(data.iloc[rr][2:],Lsize)
#         else:
#             upti = restore_upper_triangle(data.iloc[rr][2:],Lsize)
#         # print(upti)
#         if(rr%2 == 0):
#             long_res.append(upper_triangle_statics(upti,Lsize,option))
#         else:
#             short_res.append(upper_triangle_statics(upti,Lsize,option))
#     # return 1, 2
#     return merge_upti_data(long_res,option), merge_upti_data(short_res,option)


def merge_upti_data(ls,option=[]):
    res = np.zeros_like(ls[0])
    if 'err' in option:
        for data in ls:
            for i in range(data.shape[0]):
                res[i] += (data[i])**2
        return (res/len(ls))**(0.5)
    else:
        for data in ls:
            for i in range(data.shape[0]):
                res[i] += data[i]
        return res/len(ls)


def restore_1Dvec_to_matrix(data, Lsize):
    res = np.zeros([Lsize,Lsize])
    cnt = 0
    for i in range(Lsize):
        for j in range(Lsize):
            res[i][j] = data[cnt]
            cnt += 1
    return res

def restore_upper_triangle(data, Lsize):
    res = np.zeros([Lsize,Lsize])
    cnt = 0
    for i in range(Lsize):
        for j in range(Lsize-i):
            # print(i,j+i)
            res[i][j+i] = data[cnt]
            res[j+i][i] = res[i][j+i]
            cnt += 1
    return res
    # tt2 = np.zeros([6])
    # for i in range(6):
    #     for j in range(6):
    #         tt2[i] += tt[(i+j)%6][j]
    #         print(i,(i+j)%6,j)
    # print(tt)
    # print(tt2/6)

# def upper_triangle_statics(mat,Lsize,option = []):
#     if 'err' in option:
#         res = np.zeros([Lsize])
#         for i in range(Lsize):
#             for j in range(Lsize):
#                 res[i] += (mat[(i+j)%Lsize][j])**2
#         return (res/Lsize)**(0.5)
#     else:
#         res = np.zeros([Lsize])
#         for i in range(Lsize):
#             for j in range(Lsize):
#                 res[i] += mat[(i+j)%Lsize][j]
#         return res/Lsize

def upper_triangle_statics(mat,Lsize,option = []):
    if 'err' in option:
        res = np.zeros([Lsize])
        for i in range(Lsize):
            # for j in range(Lsize):
                # res[i] += (mat[(i+j)%Lsize][j])**2
            res[i] += np.std(np.diag(mat,k=i))
        return res
    else:
        res = np.zeros([Lsize])
        for i in range(Lsize):
            # for j in range(Lsize):
                # res[i] += mat[(i+j)%Lsize][j]
            res[i] += np.mean(np.diag(mat,k=i))
        return res

# def upper_triangle_statics(mat,Lsize,option = []):
#     res = np.zeros([Lsize])
#     for i in range(Lsize):
#         res[i] += mat[i][0]
#     return res

def refine_parsing_dict(parsing_dict):
    for key in parsing_dict.keys():
        all_ls = parsing_dict[key]

# def merge_csv(ls,mcs=[]):
#     res = ls[0].copy()
#     # res = pd.DataFrame(index=ls[0].index,columns=ls[0].columns)
#     # res.fillna(0)
#     # for col in res.columns:
#         # df[col].values[:] = 0

#     for col in res.columns.tolist():
#         # print(res[col])
#         res[col].values[:] = 0
#         if col.find("err") == -1:
#             for data in ls:
#                 res[col] += (data[col]/len(ls))
#         else:
#             for data in ls:
#                 res[col] += (data[col]/len(ls))/len(ls)**(0.5)
#     return res

def merge_csv(ls,mcs=[], option = False):
    res = ls[0].copy()
    # res = pd.DataFrame(index=ls[0].index,columns=ls[0].columns)
    # res.fillna(0)
    # for col in res.columns:
        # df[col].values[:] = 0

    if not mcs:
        for col in res.columns.tolist():
            # print(res[col])
            res[col].values[:] = 0
            if col.find("err") == -1:
                for data in ls:
                    res[col] += (data[col]/len(ls))
            elif option is True:
                # for data in ls:
                res[col] += ((ls[-1][col]))**2/len(ls)
                res[col] **= 0.5
            else:
                for data in ls:
                    res[col] += ((data[col])/len(ls))**2
                res[col] **= 0.5
        return res
    else :
        sum = 0
        for mm in mcs: sum += mm
        for col in res.columns.tolist():
            # print(res[col])
            res[col].values[:] = 0
            if col.find("err") == -1:
                for data, mm in zip(ls,mcs):
                    res[col] += (mm*data[col]/sum)
            else:
                for data, mm in zip(ls,mcs):
                    res[col] += (mm*(data[col])/sum)**2
                res[col] **= 0.5
        return res

def merge_dataset(parsing_dict, alpha_list):
    new_parsing_dict ={}
    for alpha in alpha_list:
        # ls, Lsize, Bins, mcs =
        a = parsing_dict[alpha]
        ls, Lsize, Bins, mcs = a[0], a[1], a[2], a[3]
        Luniq = pd.array(Lsize).unique()
        idx = 0
        res = [[],[],[],[]]
        for ll in Luniq:
            cnt = Lsize.count(ll)
            res[1].append(ll)
            tt = merge_csv(ls[idx:idx+cnt], mcs[idx:idx+cnt])
            res[0].append(tt)
            res[2].append(Bins[0])
            res[3].append(sum(mcs[idx:idx+cnt]))
            idx += cnt
        new_parsing_dict[alpha] = res
    return new_parsing_dict