"""Refresh SISIM-dependent figures using unchanged legacy combination/UQ functions.

Default: redraw from bundled small figure data (no simulation or CV).
To rebuild the data, supply --source-run, --sisim-realizations, --terrain-dir,
and --cache-dir (external generated-output directory). Requires legacy deps.
Optional --vtk writes a fresh below-ground VTK into that cache, not the repo.
"""
from pathlib import Path
import argparse, ast, contextlib, hashlib, io, json, sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'pics/final_20261018/data'
NAMES=['s_type','fs_type','p_clay','p_sand','p_gravel','entropy','var']
COLORS=['#9955bb','#f28e2b','#ffdf00'];LABELS=['Clay','Sand','Quaternary']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def extract(a):
    # Load the exact three legacy numerical functions, without importing the
    # unrelated optional GIS/fitting packages at the top of UTILITIES_Wang.
    from tqdm import tqdm
    source=ROOT/'scripts/UTILITIES_Wang.py'
    tree=ast.parse(source.read_text(encoding='utf-8-sig'))
    names={'to_3D','combine_sim_ohnelm','evaluate_nmodel'}
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==len(names)
    namespace={'np':np,'tqdm':tqdm}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),namespace)
    U=SimpleNamespace(**{k:namespace[k] for k in names})
    cache=a.cache_dir.resolve();cache.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((DATA/'exponential.json').read_text());g=cfg['preparation']['grid']
    nx,ny,nz=[g[k] for k in ['nx','ny','nz']];shape=(nz,ny,nx);ncell=nx*ny*nz
    x=g['xmn']+g['xsiz']*np.arange(nx);y=g['ymn']+g['ysiz']*np.arange(ny);z=g['zmn']+g['zsiz']*np.arange(nz)
    ix,iy=60,30;real=[0,49,99]
    q=np.load(a.source_run/'03/qbasis_realizations.npy',mmap_mode='r').reshape(100,ny,nx)
    sis=np.load(a.sisim_realizations,mmap_mode='r').reshape(100,ncell)
    terrain=np.load(a.terrain_dir/'below_gok.npy').reshape(shape);dem=np.load(a.terrain_dir/'gok.npy').reshape(ny,nx)
    dest=cache/'combined_realizations.npy'
    if dest.exists():raise FileExistsError('Use a fresh cache directory for the selected ensemble')
    model=np.lib.format.open_memmap(dest,mode='w+',dtype=np.uint8,shape=(100,ncell))
    with contextlib.redirect_stdout(io.StringIO()):
        for i in range(100):
            # Exact geometry of legacy to_3D, with its x-fastest ordering.
            qb=(q[i][None,:,:]>z[:,None,None]).ravel()
            if i in real:
                np.testing.assert_array_equal(qb,U.to_3D(q[i].ravel(),g,1)[0])
            model[i]=U.combine_sim_ohnelm(sis[i:i+1],[qb],1)[0]
    model.flush();print('Combined all 100 new SISIM fields with unchanged SGSIM surfaces',flush=True)
    out=cache/'04';out.mkdir()
    arrays={k:np.lib.format.open_memmap(out/f'{k}.npy',mode='w+',dtype=float,shape=(ncell,)) for k in NAMES}
    counts=[20,40,60,80,100];sums=np.zeros(5);n=0
    with contextlib.redirect_stdout(io.StringIO()):
        for start in range(0,ncell,20000):
            end=min(start+20000,ncell);block=model[:,start:end]
            result=U.evaluate_nmodel(block,100)
            for k,v in zip(NAMES,result):arrays[k][start:end]=v
            valid=terrain.ravel()[start:end]
            if valid.any():
                for j,count in enumerate(counts):
                    entropy=result[5] if count==100 else U.evaluate_nmodel(block[:count],count)[5]
                    sums[j]+=entropy[valid].sum()
                n+=int(valid.sum())
    for v in arrays.values():v.flush()
    frame=pd.DataFrame({'number_of_realizations':counts,'mean_entropy':sums/n,'voxel_count':n})
    frame.to_csv(DATA/'mean_entropy_by_realizations.csv',index=False)
    k=int(np.argmin(abs(z-np.median(q.mean(axis=0)))))
    bundle=dict(x=x,y=y,z=z,dem_x=dem[:,ix],dem_y=dem[iy],mask_x=terrain[:,:,ix],mask_y=terrain[:,iy,:],
        mask_plan=terrain[k],binary_sections=sis.reshape(100,*shape)[real,:,iy,:]+1,
        combined_sections=model.reshape(100,*shape)[real,:,iy,:],base_lines=q[real,iy,:])
    for name,v in arrays.items():
        field=v.reshape(shape)
        bundle[name+'_x']=field[:,:,ix];bundle[name+'_y']=field[:,iy,:];bundle[name+'_plan']=field[k]
    np.savez_compressed(DATA/'final_model_sections.npz',**bundle)
    record=dict(candidate='response_step_50',realizations=100,source_sisim_sha256=sha(a.sisim_realizations),
        source_qbasis_sha256=sha(a.source_run/'03/qbasis_realizations.npy'),combined_sha256=sha(dest),
        terrain_sha256=sha(a.terrain_dir/'below_gok.npy'),mask_count=int(terrain.sum()),
        section_x_requested=691000,section_x_center=float(x[ix]),section_y_requested=5334800,section_y_center=float(y[iy]),
        plan_elevation=float(z[k]),realization_indices=[1,50,100],
        entropy='Unchanged evaluate_nmodel, natural logs, zeros replaced by .0001, no normalization',
        classes=dict(zip(['1','2','3'],LABELS)),core_utilities_sha256=sha(ROOT/'scripts/UTILITIES_Wang.py'),
        cv='Not rerun; preserved CV figures belong to original/20261016 models, not this candidate',
        minimum_entropy_below_ground=float(arrays['entropy'][terrain.ravel()].min()),maximum_entropy_below_ground=float(arrays['entropy'][terrain.ravel()].max()))
    if a.vtk:
        import pyvista as pv
        mesh=pv.ImageData(dimensions=(nx+1,ny+1,nz+1),spacing=(g['xsiz'],g['ysiz'],g['zsiz']),
            origin=(x[0]-g['xsiz']/2,y[0]-g['ysiz']/2,z[0]-g['zsiz']/2))
        for name,v in arrays.items():mesh.cell_data[name]=np.asarray(v)
        mesh.cell_data['grid_cell_id']=np.arange(ncell)
        mesh.cell_data['DEM']=np.tile(dem.ravel(),nz)
        mesh.cell_data['entropy_class']=np.where(arrays['entropy']<.3,1,np.where(arrays['entropy']<=.7,2,3))
        mesh.cell_data['var_class']=np.where(arrays['var']<.2,1,np.where(arrays['var']<=.4,2,3))
        below=mesh.extract_cells(terrain.ravel());path=cache/'model_gok_sisim_final_20261018.vtk'
        from vtk import vtkUnstructuredGridWriter
        writer=vtkUnstructuredGridWriter();writer.SetFileName(str(path));writer.SetInputData(below)
        writer.SetFileVersion(42);writer.SetFileTypeToBinary();assert writer.Write()==1
        check=pv.read(path);ids=check.cell_data['grid_cell_id']
        assert check.n_cells==int(terrain.sum())
        for key in NAMES:np.testing.assert_array_equal(check.cell_data[key],arrays[key][ids])
        record['vtk_sha256']=sha(path);record['vtk_cells']=check.n_cells
    (DATA/'final_model_provenance.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print('New legacy UQ, plot data and optional VTK completed',flush=True)

def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(ROOT/f'pics/final_20261018/{name}.{ext}',dpi=220)
    plt.close(fig)

def plots():
    a=np.load(DATA/'final_model_sections.npz');x,y,z=[a[k] for k in ['x','y','z']]
    plt.rcParams.update({'pdf.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    cat=ListedColormap(COLORS);norm=BoundaryNorm([.5,1.5,2.5,3.5],3)
    bcat=ListedColormap(COLORS[:2]);bnorm=BoundaryNorm([.5,1.5,2.5],2)
    p=json.loads((DATA/'final_model_provenance.json').read_text())
    def image(ax,v,horizontal=x,mask=None,categorical=False,binary=False,label='',vmax=None,plan=False,tier=False):
        v=np.where(mask,v,np.nan) if mask is not None else v
        opts=dict(cmap=bcat if binary else cat,norm=bnorm if binary else norm) if categorical else dict(cmap=ListedColormap(['#70b7a4','#f5d76e','#c85a54']),norm=norm) if tier else dict(cmap='magma' if 'Entropy' in label or 'Variability' in label else 'viridis',vmin=0,vmax=vmax)
        extent=[x[0]-5,x[-1]+5,y[0]-5,y[-1]+5] if plan else [horizontal[0]-5,horizontal[-1]+5,z[0]-.25,z[-1]+.25]
        im=ax.imshow(v,origin='lower',extent=extent,aspect='equal' if plan else 'auto',interpolation='nearest',**opts)
        ax.set(xlabel='X [m]' if plan or np.array_equal(horizontal,x) else 'Y [m]',ylabel='Y [m]' if plan else 'Elevation [m]',title=label)
        ax.ticklabel_format(style='plain',useOffset=False);ax.xaxis.set_major_locator(plt.MaxNLocator(4));ax.tick_params(axis='x',labelrotation=20)
        cb=ax.figure.colorbar(im,ax=ax,shrink=.8,pad=.02)
        if categorical:cb.set_ticks([1,2] if binary else [1,2,3],labels=LABELS[:2] if binary else LABELS)
        if tier:cb.set_ticks([1,2,3],labels=['Low','Moderate','High'])
    fig,axs=plt.subplots(3,1,figsize=(11,8),layout='constrained')
    for j,ax in enumerate(axs):image(ax,a['binary_sections'][j],categorical=True,binary=True,label=f'Raw SISIM | realization {[1,50,100][j]} | Y = 5334805 m')
    save(fig,'10_sisim_sections')
    fig,axs=plt.subplots(3,1,figsize=(11,8),layout='constrained')
    axs[0].plot(x,a['base_lines'][0],label='Unchanged SGSIM Quaternary base');axs[0].plot(x,a['dem_y'],c='.3',label='Terrain')
    axs[0].set(xlabel='X [m]',ylabel='Elevation [m]',xlim=(x[0]-5,x[-1]+5));axs[0].ticklabel_format(style='plain',useOffset=False);axs[0].legend(fontsize=8)
    image(axs[1],a['binary_sections'][0],categorical=True,binary=True,label='New SISIM | realization 1')
    image(axs[2],a['combined_sections'][0],mask=a['mask_y'],categorical=True,label='New hierarchical combination | realization 1')
    axs[2].plot(x,a['dem_y'],c='.2',lw=1)
    axs[2].set_ylim(z[0]-.25,np.ceil(a['dem_y'].max()+2))
    save(fig,'11_hierarchical_combination')
    fig,axs=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    fields=[('s_type','Most likely class'),('p_clay','P(clay)'),('p_sand','P(sand)'),('p_gravel','P(Quaternary)'),('entropy','Entropy'),('var','Variability')]
    for ax,(key,label) in zip(axs.flat,fields):image(ax,a[key+'_plan'],mask=a['mask_plan'],categorical=key=='s_type',plan=True,label=label,vmax=p['maximum_entropy_below_ground'] if key=='entropy' else 1)
    fig.suptitle(f"New 100-realization ensemble | elevation {p['plan_elevation']:.2f} m")
    save(fig,'12_uq_plan')
    for key,label,lo,hi,name in [('entropy','Entropy',.3,.7,'13_sections_entropy'),('var','Variability',.2,.4,'14_sections_variability')]:
        fig,axs=plt.subplots(2,3,figsize=(16,8),layout='constrained')
        for row,(suffix,h,surface,cut) in enumerate([('x',y,a['dem_x'],'X = 691000 m'),('y',x,a['dem_y'],'Y = 5334800 m')]):
            mask=a['mask_'+suffix];v=a[key+'_'+suffix];tier=np.where(v<lo,1,np.where(v<=hi,2,3))
            image(axs[row,0],a['s_type_'+suffix],h,mask,categorical=True,label='Most likely class | '+cut)
            image(axs[row,1],v,h,mask,label=label+' | '+cut,vmax=p['maximum_entropy_below_ground'] if key=='entropy' else 2/3)
            image(axs[row,2],tier,h,mask,tier=True,label=f'{label} classes: <{lo}, {lo}-{hi}, >{hi}')
            for ax in axs[row]:ax.plot(h,surface,c='.2',lw=1);ax.set_ylim(z[0]-.25,np.ceil(surface.max()+2))
        save(fig,name)
    df=pd.read_csv(DATA/'mean_entropy_by_realizations.csv');fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    ax.plot(df.number_of_realizations,df.mean_entropy,'o-',color='#355f8d')
    ax.set(title='Mean Entropy vs. Number of Realizations',xlabel='Number of Realizations',ylabel='Mean Entropy',xticks=df.number_of_realizations,ylim=(0,df.mean_entropy.max()*1.15))
    save(fig,'16_realization_stability')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source-run','sisim-realizations','terrain-dir','cache-dir']:p.add_argument('--'+name,type=Path)
    p.add_argument('--vtk',action='store_true');a=p.parse_args()
    if any([a.source_run,a.sisim_realizations,a.terrain_dir,a.cache_dir]):
        if not all([a.source_run,a.sisim_realizations,a.terrain_dir,a.cache_dir]):p.error('Supply all four data-extraction paths')
        extract(a)
    elif a.vtk:p.error('--vtk requires the four extraction paths')
    plots()

if __name__=='__main__':main()
