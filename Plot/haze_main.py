#!/usr/bin/env python
# coding: utf-8

# 2022-02-21 v1 Try to use OOP, for Haze manuscript

# # Settings

# ## import

# In[1]:


import sys, os
import file_tool
import plot_tool
import numpy as np

from datetime import datetime
import matplotlib as mpl


# ## mpl settings

# In[2]:


fontsize = 18
mpl.rcParams.update(
    {
     'legend.fontsize': fontsize,
     'lines.linewidth' : 3,
     'axes.labelsize': fontsize+2,
     'axes.titlesize':'x-large',
     'xtick.labelsize':fontsize,
     'ytick.labelsize':fontsize,
     'font.sans-serif':'Arial',
     'font.family':'sans-serif'})


# In[ ]:





# # Function

# ## Class

# ### ?SAM

# In[3]:


class SAM:
    def __init__(self, version, case, legend, color, bin_num, style):
        self.version = version
        self.case = case
        self.style = style
        self.legend = legend
        self.bin = bin_num
        self.color = color
        
    def update(self):
        self.time_3D = [30,60]
        self.time_series = [0,60]
        self.time_z = 30
        self.path = os.path.join(r'G:\data\ACDC','v'+str(self.version),'Case'+str(self.case))
        
    def read_time_series_data(self, var_list, var_convert):
        var_new = []   
        var_calc = []
        dict_new = {}
        for var in var_list:
            if (var in var_convert.keys()):
                var_new.append(var_convert[var])
                dict_new[var_convert[var]] = var
            else:
                var_calc.append(var)
        # Read var
        data = file_tool.Lang_netCDF(self.file_STAT, 'VALUE', var_list=var_new)
        self.data = {}
        for ii in range(len(var_new)):
            self.data[dict_new[var_new[ii]]] = data[var_new[ii]]
            
        # Calcualte var
        # after minus the half day spinup, convert from day to minute
        self.data['time'] = (np.array(data['time'])-0.5) * 24 * 60 
        for var in var_calc:
            if (var == 'Rn'):
                varA = var_convert['Ra']
                varD = var_convert['Rd']
                data = file_tool.Lang_netCDF(self.file_STAT, 'VALUE', var_list=[varA,varD])
                self.data['Rn'] = data[varA] - data[varD]
                
    def plot_line(self, X, Y, ax):
        ax.plot(X, Y, color=self.color, linestyle=self.style, label=self.legend, linewidth=self.linewidth)
        
    def plot_t(self, ax, var_name):
        Y = np.mean(self.data[var_name], 1)
        self.plot_line(self.data['time'], Y, ax)
        self.calc_average(self.data['time'], Y, var_name)
        ax.set_ylim(time_lim_list[var_name])
        ax.set_ylabel(label_list[var_name])
        ax.set_xlabel(label_list['time'])
        ax.set_xlim(self.time_series)
        
    def plot_z(self, ax, var_name):
        ind_time = np.where(self.data['time']>=self.time_z)[0]
        X = np.mean(self.data[var_name][ind_time,:], 0)
        self.plot_line(X, self.data['z'], ax)
        ax.set_ylim(height_lim)
        ax.set_xlim(z_lim_list[var_name])
        ax.set_xlabel(label_list[var_name])
        ax.set_ylabel(label_list['z'])
        ax.set_xticks(z_lim_list[var_name])
        ax.set_xticklabels([str(ll) for ll in z_lim_list[var_name]])
#         ax.set_xlim(self.time_series)

    def plot_dist(self, ax):
        ax.plot(self.X_dist, self_Y_dist, )
        
    def calc_average(self, X, Y, varname):
        t1 = 30
        t2 = 60
        res = np.mean(Y[np.where((X>=t1) & (X<t2))[0]])
        print(varname,  self.legend, '30-60 ave:', res, '; 59.5:', Y[-1])


# ### ?CL

# In[4]:


class CL(SAM):
    def __init__(self, version, case, legend, color, bin_num=33, style='-'):
        legend = 'B'+legend
        super().__init__(version, case, legend, color, bin_num, style)
        self.update()
        
    def update(self):
        super().update()
        # Get File List
        step_interval = 15000
        core = '16'
        time_step = 0.02
        time_interval = int(1/time_step)
        
        self.time_list = np.arange(self.time_3D[0]*time_interval*60, self.time_3D[1]*time_interval*60+1, step_interval)
        self.file_STAT = os.path.join(self.path,'PiChamber_pi_chamber_workshop_19K.nc')
        self.file_3D = [os.path.join(self.path,'PiChamber_pi_chamber_workshop_19K_micro_'+core+'_'+str(tt).zfill(10)+'.nc') for tt in self.time_list]
        self.file_3Dc = [os.path.join(self.path,'PiChamber_pi_chamber_workshop_19K_'+core+'_'+str(tt).zfill(10)+'.nc') for tt in self.time_list]
        self.linewidth = 3
        
    def read_time_series_data(self, var_list):
        var_convert = {'T':'BTABS','Qv':'BQV','LWC':'BQC','Nd':'BNC','Na':'BNH','Reff':'BREFF',                       'Ra':'BARATE','Rd':'BDRATE','RH':'BRELH','sig':'BRSIG','time':'time','z':'z',                      'rr':'BRMEAN'}
        if (self.bin==33):
            var_convert['Na'] = 'BNA'
        super().read_time_series_data(var_list,var_convert)
     
    def read_dist_data(self):
        # X value
        var = 'RC'
        ratio = mindex**(1/3)
        temp = file_tool.Lang_netCDF(self.file_3D[0],'VALUE',var_list=[var])
        x0 = temp[var][0] # get left boundary of the 1st bin
        x1_list = [x0 * ratio**(ii)  for ii in range(self.bin)]
        x1_list = np.array(x1_list) * 10**4   # cm -> um 
        x2_list = x1_list * ratio
        bs = x2_list - x1_list
        self.X_dist = (x1_list + x2_list)/2.0        
        
        # Y value
        total_res = np.zeros((self.bin))
        var_3D_list = ['FNC'+str(ii+1).zfill(3) for ii in range(self.bin)]
        for ff in self.file_3D:
            temp = file_tool.Lang_netCDF(ff, 'VALUE', var_list=var_3D_list)
            for (ii,var) in enumerate(var_3D_list):
                total_res[ii] = total_res[ii] + np.mean(temp[var])
        N = np.array(total_res)/len(self.file_3D)
        self.Y_dist = N/bs 


# ### ?LCM

# In[5]:


class LCM(SAM):
    def __init__(self, version, case, legend, color, bin_num=33, style='-'):
        legend = 'L'+legend
        super().__init__(version, case, legend, color, bin_num, style)
        self.update()
        
    def update(self):
        super().update()
        # Get File List
        self.time_list = [self.time_3D[0]*60-1,self.time_3D[1]*60+1]
        self.file_STAT = os.path.join(self.path,'PiChamber_pi_chamber_simulation_19K_lcm.nc')
        self.file_3D = os.path.join(self.path,'PiChamber_pi_chamber_simulation_19K_lcm.spec')
        
        step_interval = 15000
        core = '16'
        time_step = 0.02
        time_interval = int(1/time_step)
        self.time_listc = np.arange(self.time_3D[0]*time_interval*60, self.time_3D[1]*time_interval*60+1, step_interval)
        self.file_3Dc = [os.path.join(self.path,'PiChamber_pi_chamber_simulation_19K_lcm_'+core+'_'+str(tt).zfill(10)+'.nc') for tt in self.time_listc]
        self.linewidth = 4
    
    def read_time_series_data(self, var_list):
        var_convert = {'T':'BTABS','Qv':'BQV','LWC':'lag_QC','Nd':'lag_NCC','Na':'lag_NAC','Reff':'lag_REF',                       'Ra':'lag_ACT','Rd':'lag_DEA','RH':'RELH','sig':'lag_RSD','time':'time','z':'z',                      'rr':'lag_R'}
        super().read_time_series_data(var_list,var_convert)
        
    def read_dist_data(self):
        ind_bin_num = 4 # line 5: how many bins
        ind_r_mean = ind_bin_num + 1 # line 6: r_mean list

        f = open(self.file_3D, 'r') 
        Lines = f.readlines() 
        LCM_bin_num = int(Lines[ind_bin_num])
        
        X = [float(xx) for xx in Lines[ind_r_mean].split()[1:]]
        self.X_dist = np.array(X) * 10**6 # m->um 

        total_res = np.zeros((LCM_bin_num))
        count = 0
        last_time = 0
        # if the first value is not time,  it would be a negative value
        for line in Lines[ind_r_mean+1:]: 
            temp = [float(jj) for jj in line.split()]
            this_time = temp[0]
            if (this_time == last_time):
                continue
            else:
                last_time = this_time
            if (this_time>self.time_list[0] and this_time<self.time_list[1]):
                total_res = total_res + temp[1:]
                count = count + 1
        res= total_res / count
        # LCM results is N/bs directly
        self.Y_dist = res / 10**12 # 1/m3/m -> 1/cm3/um


# In[6]:


#         ind_r_min = ind_bin_num + 2 # line 7: r_min list
#         ind_r_max = ind_bin_num + 3 # line 8: r_max list
#         r_min = [float(xx) for xx in Lines[ind_r_min].split()[1:]]
#         r_max = [float(xx) for xx in Lines[ind_r_max].split()[1:]]
#         bs =(np.array(r_max) - np.array(r_min))
#         bs = bs * 10**6 * 2 # m->um r->d


# ## Plots

# ### !plot_time_series

# In[7]:


def plot_time_series(obj_list, var_list, fig_ind):
    var_num = len(var_list)
    for obj in obj_list:
        obj.read_time_series_data(var_list+['time','z'])
        
    fig = plt.figure(figsize=(15,4*(var_num)))
    total_col = 6 # z plot only have 1/6 width
    for ii in range(var_num):
        axZ = plt.subplot2grid((var_num, total_col),(ii, 0))
        axT = plt.subplot2grid((var_num, total_col),(ii, 1), colspan=total_col - 1)
         
        for obj in obj_list:
            obj.plot_t(axT, var_list[ii])
            obj.plot_z(axZ, var_list[ii])
        plot_tool.label_panel(axT, label_abc[ii], pad_value=0.02, fz = 20)  
        if (ii==0):
            change_legend_order(axT,order_list)
    plt.tight_layout()
    plt.savefig(os.path.join(path_res, str(fig_ind)+'_'+date_str+'.png'),dpi=set_dpi)


# ### !plot_dist

# In[8]:


def plot_dist(obj_list, fig_ind, if_clean):
    for obj in obj_list:
        obj.read_dist_data()
        
    
    fig = plt.figure(figsize=(7,6))
    ax = plt.gca()
    for obj in obj_list:
        obj.plot_line(obj.X_dist, obj.Y_dist, ax)
        if (obj.legend=='BHsc'):
            X = obj.X_dist
            Y = obj.Y_dist
     
        
    if (dist_log[0]):
        ax.set_xscale('log')
        xlim = dist_lim_list['xlog']
        ax.set_xlim(xlim)
    else:
        xlim = dist_lim_list['xlinear']
        ax.set_xlim(xlim)
        
    if (dist_log[1]):
        ax.set_yscale('log')
        ylim = dist_lim_list['ylog']
        ax.set_ylim(ylim)
    else:
        ylim = dist_lim_list['ylinear']
        ax.set_ylim(ylim)  
        
    if (if_clean):
        ax.plot([0.9,0.9],ylim,'k:',linewidth=1)
    else:
        ax.plot([0.9,0.9],ylim,'k:',linewidth=1)
        ax.plot([0.625,0.625],ylim,':',color='red',linewidth=1)
    
        
    ax.set_xlabel(label_list['rr'])
    ax.set_ylabel(label_list['dist'])
    change_legend_order(ax,order_list)
    plt.tight_layout()
    plt.savefig(os.path.join(path_res, str(fig_ind)+'_'+date_str+'.png'),dpi=set_dpi)
    
    return X, Y


# ### !change_legend_order

# In[9]:


def change_legend_order(ax,order_list=None):
    if (order_list):
        handles, labels = ax.get_legend_handles_labels()
        nl = []
        nh = []
        for oo in order_list:
            nl.append(labels[oo])
            nh.append(handles[oo])
        ax.legend(nh, nl)
    else:
        ax.legend()


# ### obj_list

# In[10]:


def get_obj_list(if_clean=True):
    if (if_clean):
        C1 = CL(6,1000,'C','#DEB841',style='--')
        C2 = CL(6,1001,'Creg','#BF9722')
        C3 = CL(6,1002,'H','#F6A38E',40,'--')
        C4 = CL(6,1003,'Hsc','#F06543',40)
        L1 = LCM(7,1000,'','#7692FF',style='--')
        L2 = LCM(7,1001,'sc','#084887')
        # obj_list = [C1,C2,C3,C4,L1,L2]
        obj_list = [L2,L1,C2,C1,C4,C3]
    else:
        obj_list = get_obj_list(if_clean=True)
        for obj in obj_list:
            obj.case = obj.case + 10
            obj.update()
    return obj_list


# In[11]:


def get_obj_list_temp(if_clean=True):
    if (if_clean):
#         C1 = CL(6,1000,'C','#DEB841',style='--')
#         C2 = CL(6,1001,'Creg','#BF9722')
#         C3 = CL(6,1002,'H','#F6A38E',40,'--')
        C4 = CL(6,1003,'Hsc','#F06543',40)
#         L1 = LCM(7,1000,'','#7692FF',style='--')
        L2 = LCM(7,1001,'sc','#084887')
        obj_list = [L2,C4]
    else:
        obj_list = get_obj_list_temp(if_clean=True)
        for obj in obj_list:
            obj.case = obj.case + 10
            print(obj.case)
            obj.update()
    return obj_list


# ### constant

# In[12]:


def get_constant(if_clean=True):
    global label_abc, path_res, date_str, set_dpi, label_list, mindex, obj_list
    
    label_abc = ['a)','b)','c)','d)']
    # General
    path_res = r'D:\Dropbox\python\ACDC\res\GRL.v0'
    date_str = datetime.today().strftime('%m.%d_%H.%M')
    set_dpi = 300
    

    # Time Series
    label_list = {'T':'T (K)','Qv':'Q$_v$ (g/kg)', 'LWC':'LWC (g/kg)',                  'Nd':'N$_d$ (#/mg)','Na':'N$_a$ (#/mg)','Reff':'r$_{eff}$ ($\mu$m)',                  'Ra':'R$_{act}$ (#/mg/s)',                  'Rd':'R$_{dea}$ (#/mg/s)',                  'Rn':'R$_{net}$ (#/mg/s)',                  'rr':'r ($\mu$m)',                  'z':'Height (m)','time': 'Time (min)',                  'dist':'n (cm$^{-3}$$\mu$m$^{-1}$)',                  'RH':'RH (%)',                  'sig':'$\sigma_r$ ($\mu$m)',                  'diameter':'Diameter ($\mu$m)',                  'radius':'Radius ($\mu$m)'}
    mindex = 2 # double mass bin
    
    obj_list = get_obj_list(if_clean)


# # Run

# In[13]:


def run(if_clean=True):
    global time_lim_list, height_lim, z_lim_list,             order_list,              dist_log, dist_lim_list
    clean_str = "clean" if if_clean else "polluted"
    print('====================Run for ',clean_str,' cases====================')
    
    get_constant(if_clean)
    order_list = [3,2,5,4,1,0]
#     order_list = [1,0]
    


    height_lim = [0,1]

    # Distribution

    dist_log = [True, True] # if log
#     dist_log = [False,True]
    
    
    if if_clean:
        # Plot clean cases; freeze this block when plot polluted ones
        fig_diff = 0

        time_lim_list = {'T': [288,289],'Qv':[10.60,12], 'LWC': [0.02,0.1],                        'Nd':[0,55],'Na':[0,20],'Reff':[5,17.5],                        'Ra':[0,3.5],'Rd':[0,3],'Rn':[-0.1,0.7],                        'RH':[99,104],'sig':[2,5],'rr':[4,12]}
        z_lim_list = {'T': [286.5,290],'Qv':[9,13], 'LWC': [0.0,0.1],                        'Nd':[0,55],'Na':[0,20],'Reff':[7,12],                        'Ra':[0,13],'Rd':[0,5],'Rn':[-1,8],                        'RH':[99,103],'sig':[2,4],'rr':[4,8]}
        dist_lim_list = {'xlog':[0.05,100],'xlinear':[0,40],'ylog':[0.01,10**4],'ylinear':[0,30]}
    else:
        # Plot polluted cases; freeze this block when plot clean ones

        fig_diff = 4
        time_lim_list = {'T': [288,289],'Qv':[10.60,12], 'LWC': [0,0.2],                        'Nd':[40,2000],'Na':[0,2500],'Reff':[2,12],                        'Ra':[0,250],'Rd':[0,250],'Rn':[-15,15],                        'RH':[99,100],'sig':[0,5],'rr':[1,7]}
        z_lim_list = {'T': [286.5,290],'Qv':[9,13], 'LWC': [0.0,0.2],                        'Nd':[0,2500],'Na':[0,2500],'Reff':[2,8],                        'Ra':[0,1000],'Rd':[0,500],'Rn':[-100,800],                        'RH':[99,102],'sig':[0,4],'rr':[1,6]}
        dist_lim_list = {'xlog':[0.05,100],'xlinear':[0,40],'ylog':[0.01,10**6],'ylinear':[0,3000]}
    
    plot_time_series(obj_list, ['T', 'Qv', 'RH', 'LWC'], 1+fig_diff)
    plot_time_series(obj_list, ['Nd','Na','rr','sig'], 2+fig_diff)
    plot_time_series(obj_list, ['Ra','Rd','Rn'], 3+fig_diff)
    
    X, Y = plot_dist(obj_list, 4+fig_diff, if_clean)


# # main

# In[14]:


if __name__=='__main__':
    run(True)
    run(False)


# # Extra test

# ## Find dist max

# In[15]:


# plt.plot(X, Y)
# ind = 3
# print(X[ind],Y[ind])
# plt.scatter(X[ind],Y[ind],c='r')
# ax = plt.gca()
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim([10**(-3),10**2])
# np.argmax(Y)

