"""Exponential variograms and existing five-fold validation plots."""
import argparse, contextlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm

p=argparse.ArgumentParser();p.add_argument('--work-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
p.add_argument('--repo',type=Path,required=True)
p.add_argument('--part',choices=['fits','plan','cv'],default='fits');a=p.parse_args()
W=a.work_dir.resolve();R=W/'run';O=a.out.resolve();O.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.repo.resolve()/'scripts'))
import UTILITIES_Wang as U
cfg=json.loads((W/'exponential.json').read_text());g=cfg['preparation']['grid'];nx,ny,nz=g['nx'],g['ny'],g['nz']
x=g['xmn']+np.arange(nx)*g['xsiz'];y=g['ymn']+np.arange(ny)*g['ysiz'];z=g['zmn']+np.arange(nz)*g['zsiz'];xx,yy=np.meshgrid(x,y)
ext=[x[0]-5,x[-1]+5,y[0]-5,y[-1]+5];xz=[ext[0],ext[1],z[0]-.25,z[-1]+.25]
valid=np.load(W/'below_gok.npy').reshape(nz,ny,nx);gok=np.load(W/'gok.npy').reshape(ny,nx)
colors=['#9955bb','#f28e2b','#ffdf00'];cat=ListedColormap(colors);labels=['Clay','Sand','Quaternary']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):fig.savefig(O/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white');plt.close(fig)
def xyaxis(ax,section=False):
    ax.set(xlabel='X [m]',ylabel='Elevation [m]' if section else 'Y [m]');ax.ticklabel_format(style='plain',useOffset=False)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4));ax.tick_params(axis='x',labelrotation=25)
def image(ax,v,label,section=False,cmap='viridis',vmax=None,classes=False):
    kw=dict(cmap=cat,norm=BoundaryNorm([.5,1.5,2.5,3.5],3)) if classes else dict(cmap=cmap,vmin=0,vmax=vmax)
    im=ax.imshow(v,origin='lower',extent=xz if section else ext,aspect='auto' if section else 'equal',interpolation='nearest',**kw)
    xyaxis(ax,section);bar=ax.figure.colorbar(im,ax=ax,shrink=.8,pad=.03,label=label)
    if classes:bar.set_ticks([1,2,3],labels=labels)
    return im

if a.part=='fits':
    # Directional exponential fits. Simulation parameters are consolidated by
    # the unchanged legacy summary functions, recorded separately in fit_record.
    for kind,cases,name in [('continuous',[('continuous',0),('continuous',1)],'06_variogram_boundary_fit'),
                            ('categorical',[('categorical_horizontal',0),('categorical_horizontal',1),('categorical_vertical',0)],'07_variogram_indicator_fit')]:
        fig,axs=plt.subplots(1,len(cases),figsize=(5*len(cases),4.4),layout='constrained')
        record=json.loads((W/'fit_record.json').read_text());variance=record['indicator_variance']
        for j,(ax,(case,direction)) in enumerate(zip(axs,cases)):
            df=pd.read_csv(R/f'02/{case}_direction_{direction}.csv');df=df[(df.npairs>0)&(df.svargval>=0)]
            observed=df.svargval.to_numpy()*(variance if 'vertical' in case else 1)
            fit=next(f['legacy_fit'] for f in record['individual_fits'] if f['case']==case and f['direction']==direction)
            assert fit[-1]=='Exponential'
            sill,ran,nug=fit[:3]
            if 'vertical' in case:sill*=variance;nug*=variance
            h=np.linspace(0,df.avsepdist.max(),400)
            ax.scatter(df.avsepdist,observed,s=15,color='#247692',label=['E–W','N–S','Vertical'][j])
            ax.plot(h,nug+(sill-nug)*(1-np.exp(-3*h/ran)),color='#b34f43',label='Exponential')
            ax.set(xlabel='Lag distance [m]',ylabel='Semivariance',xlim=(0,h[-1]),ylim=(0,None));ax.legend(fontsize=8)
        save(fig,name)
    if a.part=='fits':sys.exit()
if a.part=='cv':
    df=pd.read_csv(R/'cv/fold_accuracy.csv');names=['clay','sand','gravel']
    cols=['fold','N','accuracy']+[n+'_indicator_acc' for n in names]
    table=[]
    for _,r in df.iterrows():table.append([str(int(r.fold)),str(int(r.N)),*[f'{r[c]:.1%}' for c in cols[2:]]])
    weighted=[np.average(df[c],weights=df.N) for c in cols[2:]]
    table.append(['All',str(int(df.N.sum())),*[f'{v:.1%}' for v in weighted]])
    fig,ax=plt.subplots(figsize=(11,3.2),layout='constrained');ax.axis('off')
    tb=ax.table(cellText=table,colLabels=['Fold','N','Overall acc.','Clay indicator acc.','Sand indicator acc.','Gravel indicator acc.'],loc='center',cellLoc='center',colWidths=[.08,.10,.16,.22,.22,.22])
    tb.auto_set_font_size(False);tb.set_fontsize(10);tb.scale(1,1.7)
    for (row,col),cell in tb.get_celld().items():cell.set_edgecolor('.8');cell.set_facecolor('#e9eef2' if row==0 else 'white')
    save(fig,'18_five_fold_accuracy_table')
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for name,color,label in zip(names,colors,labels):ax.plot(df.fold,df[name+'_indicator_acc'],'o-',color=color,label=label)
    ax.set(xlabel='Fold',ylabel='Indicator accuracy',xticks=range(1,6),ylim=(0,1));ax.legend(fontsize=9);save(fig,'18_indicator_prediction_accuracy')
    df.to_csv(O/'five_fold_accuracy.csv',index=False)

if a.part in ('plan','cv'):
    s=pd.read_csv(R/'cv/borehole_folds.csv');s=s[s.x.between(ext[0],ext[1])&s.y.between(ext[2],ext[3])]
    fig,axs=plt.subplots(1,5,figsize=(17,3.8),layout='constrained')
    for i,ax in enumerate(axs):
        train=s[s.fold!=i];test=s[s.fold==i]
        ax.scatter(train.x,train.y,s=4,color='.78',label='Train');ax.scatter(test.x,test.y,s=8,color='#b34f43',label='Test')
        ax.set(xlim=ext[:2],ylim=ext[2:],title=f'Fold {i+1}');ax.set_aspect('equal');xyaxis(ax);ax.xaxis.set_major_locator(plt.MaxNLocator(2))
        if i==0:ax.legend(fontsize=8)
    save(fig,'19_five_fold_plan')
    print('CV PNGs complete',flush=True)
