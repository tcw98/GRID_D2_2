"""Redraw the final input/output variograms from the small bundled plot data.

No simulation, fitting, GAMV rerun, or cross-validation is performed here.
Run: python scripts/plot_final_variograms.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'pics/final_20261018/data'

def theory(h,model,direction):
    # GSLIB zero-angle major axis is NS; minor axis is EW.
    index={'NS':5,'EW':6,'Z':7}[direction]
    assert model['c0']==0 and all(row[0]==2 for row in model['nest'])
    return sum(row[1]*(1-np.exp(-3*np.asarray(h)/row[index])) for row in model['nest'])

def main():
    cfg=json.loads((DATA/'exponential.json').read_text())
    model=cfg['variograms']['categorical']
    norm=json.loads((DATA/'variogram_normalization.json').read_text())
    marg=pd.read_csv(DATA/'marginals.csv')
    plt.rcParams.update({'pdf.fonttype':42,'font.size':10})
    for prefix,scope in [('neogene','Actual Neogene region'),('full','Complete SISIM grid')]:
        for mode in ['shape','common_variance']:
            fig,axes=plt.subplots(1,3,figsize=(13.3,4.3),layout='constrained')
            for d,ax in zip(['EW','NS','Z'],axes):
                t=pd.read_csv(DATA/f'experimental_{d}.csv')
                t=t[(t.pairs>0)&(t.distance>0)&(t.nominal_distance>0)]
                values=np.load(DATA/f'{prefix}_{d}_curves.npy')
                assert values.shape[0]==100
                h=np.arange(1,values.shape[1]+1)*(.5 if d=='Z' else 10)
                ediv=tdiv=rdiv=norm['common_reference_variance']
                if mode=='shape':
                    p=marg['p_neogene' if prefix=='neogene' else 'p_full'].to_numpy()
                    ediv=norm['experimental_variance'];tdiv=norm['input_total_sill']
                    rdiv=(p*(1-p))[:,None]
                a=values/rdiv;lo,hi=np.quantile(a,[.05,.95],axis=0)
                hp=np.r_[0,h];limit=15 if d=='Z' else 155
                ax.scatter(t.distance,t.gamma/ediv,c='#222222',s=23,zorder=4,label='Experimental')
                ax.fill_between(hp,np.r_[0,lo],np.r_[0,hi],color='#007f86',alpha=.2,edgecolor='none',linewidth=0,label='5th-95th percentile')
                ax.plot(hp,np.r_[0,a.mean(axis=0)],'-o',c='#007f86',lw=1.8,ms=2.6,label='Mean reproduced')
                x=np.linspace(0,limit,500)
                ax.plot(x,theory(x,model,d)/tdiv,'--',c='#ce8a1a',lw=1.6,label='Calibrated theoretical input')
                ax.set(xlim=(0,limit),ylim=(0,1.2),title=d,xlabel='Separation distance [m]',ylabel='Normalized semivariance')
                ax.grid(alpha=.16)
                assert len(ax.lines)==2 and ax.lines[0].get_xdata()[0]==0 and ax.lines[0].get_ydata()[0]==0
            handles,labels=axes[0].get_legend_handles_labels(); lookup=dict(zip(labels,handles))
            order=['Experimental','Calibrated theoretical input','Mean reproduced','5th-95th percentile']
            axes[0].legend([lookup[k] for k in order],order,fontsize=7.4,loc='lower right')
            text='Reference sill = 1: separate variance/sill normalization (shape only)' if mode=='shape' else 'Common reference variance 0.24815'
            fig.suptitle(f'{scope} | 100 realizations\n{text}',fontsize=11)
            name=f'08_variogram_{prefix}_{mode}'
            for ext in ['png','pdf']:fig.savefig(ROOT/f'pics/final_20261018/{name}.{ext}',dpi=210)
            plt.close(fig)
    # Input-model panel uses the same target bins and normalization as the shape plot.
    fig,axes=plt.subplots(1,3,figsize=(13.3,4.3),layout='constrained')
    for d,ax in zip(['EW','NS','Z'],axes):
        t=pd.read_csv(DATA/f'experimental_{d}.csv');t=t[(t.pairs>0)&(t.distance>0)&(t.nominal_distance>0)]
        limit=15 if d=='Z' else 155;x=np.linspace(0,limit,500)
        ax.scatter(t.distance,t.gamma/norm['experimental_variance'],c='#222222',s=23,label='Experimental')
        ax.plot(x,theory(x,model,d)/norm['input_total_sill'],'--',c='#ce8a1a',label='Calibrated theoretical input')
        ax.set(xlim=(0,limit),ylim=(0,1.2),title=d,xlabel='Separation distance [m]',ylabel='Normalized semivariance')
        ax.grid(alpha=.16)
    axes[0].legend(fontsize=8,loc='lower right');fig.suptitle('Final four-structure SISIM input | zero nugget | reference sill = 1')
    for ext in ['png','pdf']:fig.savefig(ROOT/f'pics/final_20261018/07_variogram_indicator_fit.{ext}',dpi=210)
    plt.close(fig)

if __name__=='__main__':main()
