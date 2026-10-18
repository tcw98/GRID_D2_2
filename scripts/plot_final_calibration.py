"""Concise Wiegel-style final SISIM and unchanged SGSIM calibration figures.

Default: redraw the three panels from bundled tables; no simulation or refit.
SISIM: 2% probability-bin staircase, raw binary held-out Neogene predictions.
SGSIM: ECDF of saved out-of-fold local quantile ranks (mid-PIT).
Neither result is smoothed. Empty binary bins remain missing.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'pics/final_20261018/data'
TEAL='#007f86';ORANGE='#d28c23'

def panel(ax,kind):
    ax.plot([0,1],[0,1],'--',c=ORANGE,lw=1.5,label='Bi-sector line',zorder=2)
    if kind=='SISIM':
        t=pd.read_csv(DATA/'binary_sisim_calibration_2pct.csv')
        assert len(t)==50
        assert t['count'].sum()==22313
        edges=np.r_[t.lower.iloc[0],t.upper.to_numpy()]
        ax.stairs(t.observed_frequency.to_numpy(),edges,baseline=None,color=TEAL,lw=1.5,label='Cross-validation',zorder=3)
        ax.set(xlabel='Predicted P(Sand)',ylabel='Observed Sand frequency',title='Binary SISIM')
    else:
        t=pd.read_csv(DATA/'continuous_heldout_predictions.csv')
        pit=np.sort(t.pit.to_numpy());n=len(pit)
        # Exact empirical CDF, including finite-ensemble rank ties.
        ax.step(np.r_[0,pit,1],np.r_[0,np.arange(1,n+1)/n,1],where='post',color=TEAL,lw=1.5,label='Cross-validation',zorder=3)
        ax.set(xlabel='Predicted quantile',ylabel='Cumulative validation frequency',title='SGSIM')
    ax.set(xlim=(0,1),ylim=(0,1),aspect='equal',xticks=np.arange(0,1.01,.2),yticks=np.arange(0,1.01,.2))
    ax.grid(alpha=.14);ax.legend(loc='upper left',frameon=False,fontsize=8)
    assert len(ax.lines)==(1 if kind=='SISIM' else 2)

def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(ROOT/f'pics/final_20261018/{name}.{ext}',dpi=240)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sgsim-only',action='store_true',help='Redraw unchanged SGSIM data without waiting for SISIM')
    a=p.parse_args()
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    if a.sgsim_only:
        fig,ax=plt.subplots(figsize=(5.1,5.0),layout='constrained');panel(ax,'SGSIM')
        fig.suptitle('Five-fold cross-validation | 100 realizations per fold',fontsize=10)
        save(fig,'figure9_sgsim_calibration');return
    metadata=json.loads((DATA/'final_cv_summary.json').read_text(encoding='utf-8'))
    assert metadata['status']=='complete' and metadata['candidate']=='response_step_50'
    assert metadata['realizations_per_fold']==[100]*5
    fig,axes=plt.subplots(1,2,figsize=(10.0,4.9),layout='constrained')
    for ax,kind in zip(axes,['SISIM','SGSIM']):panel(ax,kind)
    fig.suptitle('Five-fold cross-validation | 100 realizations per fold',fontsize=11)
    save(fig,'figure9_final_calibration')
    for kind in ['SISIM','SGSIM']:
        fig,ax=plt.subplots(figsize=(5.1,5.0),layout='constrained');panel(ax,kind)
        fig.suptitle('Five-fold cross-validation | 100 realizations per fold',fontsize=10)
        save(fig,'figure9_'+kind.lower()+'_calibration')

if __name__=='__main__':main()
