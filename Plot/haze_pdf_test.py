#!/usr/bin/env python
# coding: utf-8

# In[1]:


import phys_tool as pt
import file_tool as ft
import plot_tool as plot
import import_ipynb
import haze_main as haze


# In[2]:


def get_SAM_P(file):
    A = ft.Lang_netCDF(file,'VALUE',['p','PP'])
    shp = np.shape(A['PP'])
    new_p = np.tile(A['p'][:, np.newaxis, np.newaxis], (1, shp[1], shp[2]))
    return new_p*100 + A['PP']


# In[3]:


def find_bin_value(value, x_list):
    x0 = x_list[0]
    Y = []
    Xm = []
    for x1 in x_list[1:]:
        res = len(np.where((value>=x0) & (value<x1))[0])
#         print(x0,x1,res)
        Y.append(res)
        Xm.append((x0+x1)/2)
        x0 = x1
    return Xm, Y



# In[4]:


def bin_all_file(obj):
    RHall = []
    for file in obj.file_3Dc:
        A = ft.Lang_netCDF(file,'VALUE',[tvar,qvar])
        P = get_SAM_P(file)
        A['s'] = []
        for qq, tt, pp in zip(A[qvar].flatten(), A[tvar].flatten(), P.flatten()):
    #         print(qq,tt)
            A['s'].append(pt.supersaturation(qq/1000, tt, kind = 'SAM', P=pp))
        A['RH'] = (np.array(A['s']) + 1) * 100
#         RHall.extend(A['RH'][(A['RH']<=101)&(A['RH']>=99)])
        RHall.extend(A['RH'])
    
    temp1 = np.mean(A['RH'][A['RH']<=100.08])
    temp2 = np.mean(A['RH'][A['RH']<=100])
    print(obj.legend,'[100.08]',temp1,'[100]',temp2)
        
#         
#     x1 = np.floor(np.min(A['RH'])*10)/10
#     x2 = np.ceil(np.max(A['RH'])*10)/10
    x1 = 99
    x2 = 101
    dx = 0.01
    X = np.arange(x1,x2,dx)
    Xm, Y = find_bin_value(RHall, X)
    return Xm, Y


# In[5]:


def plot_bin_XY(obj, ax, max_Y=0):
    X, Y = bin_all_file(obj) 
    Y = np.array(Y)/1000 # Population (x1000)
    obj.plot_line(X, Y, ax)
    max_Y = np.max([max_Y, np.max(Y)])
    ind = np.argmax(Y)
    print(obj.legend, X[ind], Y[ind])
    return X, Y, max_Y


# In[6]:



qvar = 'QV'
tvar = 'TABS'
RHvar = 'BRELH'


# In[7]:


fig = plt.figure(figsize=(14,6))
for cc, ic in enumerate([True, False]):
    haze.get_constant(if_clean=ic)
    if ic:
        ax = plt.subplot(1, 2, cc+1)
    else:
        ax = plt.subplot(1, 2, cc+1)
    
#     ax = plt.gca()
    max_Y = 0

    for obj in [haze.obj_list[ii] for ii in [0,4]]:
        X, Y, max_Y = plot_bin_XY(obj, ax, max_Y)

    ax.plot([100.08,100.08], [0,max_Y*1.1], 'k:', linewidth = 1) #, label='Critical RH'
    if cc==1:
        ax.plot([100.045,100.045], [0,max_Y*1.1], 'r:', linewidth = 1) #peak
    if cc==0:
        ax.set_title('0.25 cm$^{-3}$s$^{-1}$ injection',fontsize=haze.fontsize)
    else:
        ax.set_title('2.5 cm$^{-3}$s$^{-1}$ injection',fontsize=haze.fontsize)
    ax.set_xlim([99,101])
    ax.set_xticks([99,99.5,100,100.5,101])
    ax.set_ylim([0,max_Y*1.1])
    ax.set_xlabel('RH (%)')
    ax.set_ylabel('Population (x10$^3$)')
    
    plot.label_panel(ax, haze.label_abc[cc], pad_value=0.02, fz = 20)  
    if ic:
        plt.legend()
        
plt.savefig(haze.path_res+r'/HR_pdf.png',dpi = haze.set_dpi)




