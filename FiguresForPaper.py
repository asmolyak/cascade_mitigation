


'''
Figure 1:
1. Run bipartite_cascade_analytical_solution_sergey.py
Get agg_all.
Change to pdagg = pd.DataFrame(agg_all, index=threshold)
Open results2018-08-25.pgz in bipartite_pickle_postprocess.py
run er_health, er_iters, er_avdeg, er_max_clusts = sum_up_new(a[(True,'er')])
plot as per below:
    
plt.figure(figsize=(16,12))

ax = pdagg[8].iloc[:12].plot(style='.-',linewidth=10)
(a[(True,'er')][0][0][8]/10000).plot(style = 'o',markeredgewidth=10,markersize=20,markerfacecolor='none',ax=ax)
(er_health[8].iloc[:12]/10000).plot(style='--',ax=ax)
ax = pdagg[8].iloc[12:].plot(style='.-',color='#1f77b4',linewidth=10)
(er_health[8].iloc[12:]/10000).plot(style='--',color='#2ca02c',ax=ax)
ax.legend(['Analytic Solution','Typical Realization','Multiple Iterations'],prop={'size':30})
ax.tick_params(labelsize=30)
ax.set_xlabel('Threshold',fontsize=30)
ax.set_ylabel('Surviving Fraction of the Network',fontsize=30)


fig 1 b
#fig2 = {}
#for c in pdagg.columns:
#    print c,1/(1-threshold[np.argmin(np.diff(pdagg[c].values))+1])
#ax=pd.DataFrame(fig2).T[0].plot()
fig2 = {}
fig3 = {}
for c in pdagg.columns:
    fig3[c]=threshold[np.argmin(np.diff(er_health[c].values))]
    fig2[c]=threshold[np.argmin(np.diff(pdagg[c].values))]
plt.figure(figsize=(16,12))
ax = pd.Series(fig2).plot(style='o',markersize=20,markerfacecolor='none',markeredgewidth=6)
pd.Series(fig3).plot(ax=ax,style='.',markersize=20)
ax.tick_params(labelsize=30)
ax.legend(['Analytic Solution','Simulation'],fontsize=30)
ax.set_title('Critical threshold by Initial Degree',fontsize=30)
ax.set_xlabel('Initial Degree',fontsize=30)
ax.set_ylabel('Critical Threshold',fontsize=30)

fig 1 c number of iterations per threshold

plt.figure(figsize=(16,12))
ax = er_iters[8].plot(linewidth=10)
ax.set_xlabel('Threshold',fontsize=30)
ax.set_ylabel('Average Number of Iterations',fontsize=30)
ax.tick_params(labelsize=30)


fig 1 d - iterations at criticality
max_iters = {}
for i, j in er_health.diff().fillna(0).abs().idxmax().items():
    max_iters[i] = er_iters[i][j]

plt.figure(figsize=(16,12))
ax = pd.Series(max_iters).plot(style='.',markersize=20)
ax.set_xlabel('Initial Average Degree',fontsize=30)
ax.set_ylabel('Number of Iterations',fontsize=30)
ax.tick_params(labelsize=30)

fig 1 e
er_comp_sz.columns = list(product(np.arange(5,20),1-1.0/np.arange(4.0, 17.0, 0.25)))

plt.figure(figsize=(16,12))
ax = (er_comp_sz[(8,1-1/7)][1:12]/10000).plot(linewidth=10)
ax.set_xlabel('Cascade Iteration',fontsize=30)
ax.set_ylabel('Largest Surviving Component',fontsize=30)
ax.tick_params(labelsize=30)

fig 1 f
er_crsh_sz.columns = list(product(np.arange(5,20),1-1.0/np.arange(4.0, 17.0, 0.25)))
plt.figure(figsize=(16,12))
ax = (er_crsh_sz[(8,1-1/7)][1:12]/10000).plot(linewidth=10)
ax.set_xlabel('Cascade Iteration',fontsize=30)
ax.set_ylabel('Fraction of Failure per Iteration',fontsize=30)
ax.tick_params(labelsize=30)


figure 3,4 a,c,e:
run cent_temp_correct_cascade_multi with: 
    3a: g = nx.fast_gnp_random_graph(10000,avdeg/10000), protect=['model'] and dontsave=np.arange(0,1.1,0.1)
    3c: g = sf_graph(10000,power=2.5,avdeg=avdeg), protect=['model'] and dontsave=np.arange(0,1.1,0.1)
    3e: round_rvs = generate_distribution(p=2.5,min_holdings=1,num_of_banks=5000,\
                                          num_of_assets=5000,fit=False,fixed_mean=True,average=avdeg)
        single_net = False
        project_network = False
        g = setup_network(num_of_banks=5000,num_of_assets=5000,\
                          round_rvs=round_rvs,p=2.5,kind='man_half_half',min_holdings=1,\
                          single_net=single_net,project_network=project_network,av_deg=avdeg)

                , protect=['model'] and dontsave=np.arange(0,1.1,0.1)
    4a,c,e same network structures as 3, protect = ['model','random','top_degree','bottom_deg'] and dontsave = [0] 
plot:
    fig = plt.figure(6,figsize=(16,12))
    plt.tick_params(labelsize=20)

#    plt.plot(df_final.index, df_final[df_final.columns[0]],'-.',label='0.0')
#    plt.plot(df_final.index, df_final[df_final.columns[3]],'-.',label='0.3')
#    plt.plot(df_final.index, df_final[df_final.columns[6]],'-.',label='0.6')
#    plt.plot(df_final.index, df_final[df_final.columns[9]],'-.',label='1.0')
    
    plt.errorbar(df_final.index, df_final[df_final.columns[0]], yerr=df_final[[df_final.columns[1], df_final.columns[2]]].abs().T.values,linestyle='-.',capsize=3,label='0')
    plt.errorbar(df_final.index, df_final[df_final.columns[3]], yerr=df_final[[df_final.columns[4], df_final.columns[5]]].abs().T.values,linestyle='-.',capsize=3,label='0.3')
    plt.errorbar(df_final.index, df_final[df_final.columns[6]], yerr=df_final[[df_final.columns[7], df_final.columns[8]]].abs().T.values,linestyle='-.',capsize=3,label='0.6')
    plt.errorbar(df_final.index, df_final[df_final.columns[9]], yerr=df_final[[df_final.columns[10], df_final.columns[1]]].abs().T.values,linestyle='-.',capsize=3,label='1.0')

    plt.legend(loc=6,fontsize=20)
    plt.xlabel('Threshold',fontsize=20)
    plt.ylabel('Probability of System survival',fontsize=20)

    axes1 = fig.add_subplot(111)
    axes2 = axes1.twinx()
    df_count['model_0_mean'].plot(ax=axes2,style='--r',ylim=[0.0,1.0],label='Protected Nodes')
    plt.ylabel('Fraction of nodes held safe',fontsize=20)
    plt.legend(loc=0,fontsize=20)
    plt.tick_params(labelsize=20)
    plt.show()

figure 3b:
    fig = plt.figure(6,figsize=(16,12))
    plt.tick_params(labelsize=20)
    objs = plt.plot(np.arange(0,1.1,0.1),df_final.loc[[0.865,0.885,0.895]][['model_{}_mean'.format(np.round(i,2)) for i in np.arange(0,1.1,0.1)]].T)
    plt.xlabel('Probability of a Safe node to Fail',fontsize=20)
    plt.ylabel('Probability of System survival',fontsize=20)
    plt.legend(iter(objs), ('0.865','0.885','0.895'),loc=0,fontsize=20)
    plt.tick_params(labelsize=20)

Figure 5 for paper - data visualisation
in real_data_irena.py run from line 88 onward
def get_cmap(n, name='tab20'):
#    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#    RGB color; the keyword argument name must be a standard mpl colormap name.
    return plt.cm.get_cmap(name, n)

doms = {j[:2] for j in bottom}
clrs = {}
for n,d in enumerate(doms):
    clrs[d] = n
cmap = get_cmap(len(doms))
x=[];y=[];size=[];color=[]
for i in bottom:
    x.append(sovereign_debt_copy.node[i]['init_value'])
    y.append(np.random.rand())
    size.append(sovereign_debt_copy.degree(i))
    color.append(cmap(clrs[i[:2]]))
# fig 5a 
fig = plt.figure(figsize=(16,12))
plt.scatter(x,y,s=[(kr*50) for kr in size],c=color)
plt.legend([mpatches.Patch(color=cmap(b)) for b in clrs.values()],
           [i for i in clrs.keys()])
plt.xlabel('Overall Exposure',fontsize=20)
plt.yticks([])


# fig 5b 
fig = plt.figure(figsize=(16,12))
plt.tick_params(labelsize=20)
plt.scatter(x,y,s=[(kr*50) for kr in size],c=color)
plt.semilogx(True)
plt.legend([mpatches.Patch(color=cmap(b)) for b in clrs.values()],
           [i for i in clrs.keys()])
plt.xlabel('Overall Exposure (Log scale)',fontsize=20)
plt.yticks([])

#fig 5c number of banks per country

count = {}; exposure = {}
for i in bottom:
    if i[:2] in count.keys():
        count[i[:2]] += 1
        exposure[i[:2]] += sovereign_debt_copy.node[i]['init_value']
    else:
        count[i[:2]] = 1
        exposure[i[:2]] = sovereign_debt_copy.node[i]['init_value']

fig = plt.figure(figsize=(16,12))
plt.tick_params(labelsize=20)
sr = pd.Series(count)
sr.sort_values(ascending=False,inplace=True)
sr.plot(kind='bar',color=[cmap(clrs[k[:2]]) for k in sr.index])
plt.ylabel('Number of Banks',fontsize=20)
plt.xlabel('Country',fontsize=20)


fig = plt.figure(figsize=(16,12))
plt.tick_params(labelsize=20)
exp = pd.Series(exposure)
exp.sort_values(ascending=False,inplace=True)
exp.plot(kind='bar',color=[cmap(clrs[k[:2]]) for k in exp.index],logy=True)
plt.ylabel('Exposure',fontsize=20)
plt.xlabel('Country',fontsize=20)



end fig 5
supplementary

load data:
    
with gzip.open('results2018-12-12_1000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health1, er_iters1, er_avdeg1, er_max_clusts1, er_crsh_sz1, er_comp_sz1 = sum_up_new(a[True])

with gzip.open('results2018-12-12_2000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health2, er_iters2, er_avdeg2, er_max_clusts2, er_crsh_sz2, er_comp_sz2 = sum_up_new(a[True])

with gzip.open('results2018-12-12_4000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health4, er_iters4, er_avdeg4, er_max_clusts4, er_crsh_sz4, er_comp_sz4 = sum_up_new(a[True])

with gzip.open('results2018-12-12_8000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health8, er_iters8, er_avdeg8, er_max_clusts8, er_crsh_sz8, er_comp_sz8 = sum_up_new(a[True])

with gzip.open('results2018-12-12_16000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health16, er_iters16, er_avdeg16, er_max_clusts16, er_crsh_sz16, er_comp_sz16 = sum_up_new(a[True])

with gzip.open('results2018-12-13_32000Parallel.pgz','r') as h:
    a = pickle.load(h)  # 2018-08-22, (a[(True,'er')][56][0][8]/4000).plot()
er_health32, er_iters32, er_avdeg32, er_max_clusts32, er_crsh_sz32, er_comp_sz32 = sum_up_new(a[True])


supplementary unscaled iterations:
plt.figure(figsize=(16,12))
er_iters1[8].plot()
er_iters2[8].plot()
er_iters4[8].plot()
er_iters8[8].plot()
er_iters16[8].plot()
er_iters32[8].plot()
plt.xlabel('Threshold',fontsize=30)
plt.ylabel('Number of Iterations',fontsize=30)
plt.tick_params(labelsize=30)
plt.legend(['1000 nodes','2000 nodes','4000 nodes','8000 nodes','16000 nodes','32000 nodes'],fontsize=30)

supplementary scaled iterations:
ratios = [1,2,4,8,16,32]

plt.figure(figsize=(16,12))
(((1.0/ratios[0])**alpha)*er_iters1[8]).plot()
(((1.0/ratios[1])**alpha)*er_iters2[8]).plot()
(((1.0/ratios[2])**alpha)*er_iters4[8]).plot()
(((1.0/ratios[3])**alpha)*er_iters8[8]).plot()
(((1.0/ratios[4])**alpha)*er_iters16[8]).plot()
(((1.0/ratios[5])**alpha)*er_iters32[8]).plot()
plt.legend(['1000 nodes','2000 nodes','4000 nodes','8000 nodes','16000 nodes','32000 nodes'],fontsize=30)
plt.xlabel('Threshold',fontsize=30)
plt.ylabel('Number of Iterations',fontsize=30)
plt.tick_params(labelsize=30)

supplementary scaling
run pickle postprocess to calculate max_1,2,3
linregress(np.log([1000,2000,4000,8000,16000,32000]),np.log(max_2))
line = lambda x: x*0.1326325475290079 + 0.7894845636178769
plt.figure(figsize=(16,12))
plt.plot(np.log([1000,2000,4000,8000,16000,32000]),np.log(max_2),'.',markersize=40)
plt.plot(np.log([1000,2000,4000,8000,16000,32000]),line(np.log([1000,2000,4000,8000,16000,32000])),linewidth=10)
plt.legend(['Empirical Measurements','Linear Fit, Slope=0.133'],fontsize=30)
plt.xlabel('Log Network Size',fontsize=30)
plt.ylabel('Log Number of Iterations at criticality',fontsize=30)
plt.tick_params(labelsize=30)
'''