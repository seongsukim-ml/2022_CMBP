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

alpha_file_name_dict = {
                    "a=1"    :"1.000000",
                    "a=1.5"  :"1.500000",
                    "a=2"    :"2.000000",
                    "a=2.5"  :"2.500000",
                    "a=3"    :"3.000000",
                    "a=3.5"  :"3.500000",
                    "a=100"  :"100.000000",
}

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
        Sus = (a.loc[:,'mm**2'].values) - MM**2
        MMerror = a.loc[:,'MMerr'].values
        CCerror = a.loc[:,'CCerr'].values
        # MMerror = np.sqrt(1/(18000-1)*abs(a.iloc[:,5].values/10000-(a.iloc[:,4].values/100)**2))
        Suserr  = abs(a.loc[:,'MM2err'].values-abs(MM*2*MMerror))
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

def draw_binder_FFS(path_list, Lsize, Bins, Step_Size, xliml = None, yliml = (0.75,1), Tc=2.269, nu=1):
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
    # plt.xlabel('$L^{1/\nu}[T-T_c]$')
    plt.show()

def draw_binder_multi(path_list,Lsize,Bins,Step_Size, Tc=2.269, nu=1 , plot_list = ["MM","CC"],xliml=[],yliml = []):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    marker = ["o","s","o","s","o"]

    for path in path_list:
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
            BBerr = a.loc[:,"BBerr"].values
            plt.errorbar(TT,Binder,BBerr,marker=marker[0],mfc='none',markersize=5,capsize= 7)
        else :
            plt.plot(TT,Binder,marker=marker[0],mfc='none',markersize=5)
        # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])
    if xliml:
        plt.xlim(xliml)
    if yliml:
        plt.ylim(yliml)
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.ylabel('Binder ratio g')
    plt.xlabel('Temperature')

    plt.show()

def draw_binder_FFS_multi(path_list,Lsize,Bins,Step_Size, Tc=2.269, nu=1 , plot_list = ["MM","CC"]):
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
            print("calc")
            Binder = 0.5*(3-a.loc[:,"mm**4"].values/(a.loc[:,"mm**2"].values)**2)
        else:
            Binder = a.loc[:,'Binder'].values
        # plt.ylim(-0.1,2)
        # plt.xlim(1.5,4)
        plt.style.use('seaborn-whitegrid')
        TTa = TT-Tc
        
        if "BBerr" in plot_list:
            BBerr = a.loc[:,"BBerr"].values
            plt.errorbar(TTa*(Lsize[i]**(1/nu)),Binder,BBerr,linestyle="",marker=marker[0],markersize=5,capsize=7,mfc='none',label="L"+str(Lsize[i]))
        else :
            plt.plot(TTa*(Lsize[i]**(1/nu)),Binder,linestyle="",marker=marker[0],markersize=5,mfc='none',label="L"+str(Lsize[i]))
        # plt.axvline(x=2/np.log(1+np.sqrt(2)),c='grey',lw=1,dashes=[2,2])

        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
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
    parsing_dict_cor = {}
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
        if(ls.find('cor') == -1):
            parsing_dict[dir]= [ls, Lsize, Bins, mcs]
        else:
            parsing_dict_cor[dir]= [ls, Lsize, Bins, mcs]
    return parsing_dict, parsing_dict_cor, alpha_list

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