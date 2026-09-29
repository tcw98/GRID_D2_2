"""Plot existing realizations; preserve legacy uncertainty and unrotated coordinates."""
import argparse, contextlib, io, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

p=argparse.ArgumentParser()
for name in ['work-dir','repo','out']: p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--part',choices=['stability','map','sections','vtk'],required=True)
p.add_argument('--x-cut',type=float); p.add_argument('--y-cut',type=float)
a=p.parse_args(); W=a.work_dir; R=W/'run'; O=a.out; O.mkdir(parents=True,exist_ok=True)
(O/'data').mkdir(exist_ok=True)
g=json.loads((W/'exponential.json').read_text(encoding='utf8'))['preparation']['grid']
nx,ny,nz=[g[k] for k in ['nx','ny','nz']]
x=g['xmn']+np.arange(nx)*g['xsiz']; y=g['ymn']+np.arange(ny)*g['ysiz']; z=g['zmn']+np.arange(nz)*g['zsiz']
extent=[x[0]-5,x[-1]+5,y[0]-5,y[-1]+5]
mask=np.load(W/'below_gok.npy').reshape(nz,ny,nx); dem=np.load(W/'gok.npy').reshape(ny,nx)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    fig.savefig(O/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white'); plt.close(fig)
def axis(ax,horizontal='X [m]',vertical='Y [m]'):
    ax.set(xlabel=horizontal,ylabel=vertical); ax.ticklabel_format(style='plain',useOffset=False)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4)); ax.tick_params(axis='x',labelrotation=20)
def classes(v,lo,hi): return np.where(v<lo,1,np.where(v<=hi,2,3)).astype(np.uint8)

if a.part=='stability':
    sys.path.insert(0,str(a.repo/'scripts')); import UTILITIES_Wang as U
    model=np.load(R/'03/combined_realizations.npy',mmap_mode='r')
    counts=sorted(set(range(20,len(model)+1,20))|{len(model)})
    totals=np.zeros(len(counts)); n=0
    reference=np.load(R/'04/entropy.npy',mmap_mode='r')
    with contextlib.redirect_stdout(io.StringIO()):
        for start in range(0,mask.size,20000):
            ids=np.flatnonzero(mask.ravel()[start:start+20000])+start
            if not len(ids): continue
            ent=[U.evaluate_nmodel(model[:k,ids],k)[5] for k in counts]
            np.testing.assert_allclose(ent[-1],reference[ids],rtol=0,atol=0)
            totals += [v.sum() for v in ent]; n+=len(ids)
    df=pd.DataFrame({'number_of_realizations':counts,'mean_entropy':totals/n,'voxel_count':n})
    df.to_csv(O/'data/mean_entropy_by_realizations.csv',index=False)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained'); ax.plot(counts,totals/n,'o-',color='#355f8d')
    ax.set(title='Mean Entropy vs. Number of Realizations',xlabel='Number of Realizations',ylabel='Mean Entropy',ylim=(0,float((totals/n).max())*1.15))
    if len(counts)<=15: ax.set_xticks(counts)
    save(fig,'16_realization_stability'); print(df.to_string(index=False))

if a.part=='map':
    import requests
    from PIL import Image
    from pyproj import Transformer
    from rasterio.transform import from_bounds
    from rasterio.warp import reproject, Resampling
    # Fetch only this map viewport, one zoom level, and retain the local cache.
    cache=O.parent/'osm_cache'; cache.mkdir(exist_ok=True); zoom=16; world=20037508.342789244
    to_web=Transformer.from_crs(25832,3857,always_xy=True)
    corners=np.array([to_web.transform(i,j) for i in extent[:2] for j in extent[2:]])
    xmin,ymin=corners.min(axis=0); xmax,ymax=corners.max(axis=0); span=2*world/2**zoom
    tx0,tx1=np.floor((np.array([xmin,xmax])+world)/span).astype(int)
    ty0,ty1=np.floor((world-np.array([ymax,ymin]))/span).astype(int)
    mosaic=Image.new('RGB',((tx1-tx0+1)*256,(ty1-ty0+1)*256))
    session=requests.Session(); session.headers['User-Agent']='GRID-D22-ScientificFigure/1.0'
    urls=[]
    for tx in range(tx0,tx1+1):
        for ty in range(ty0,ty1+1):
            file=cache/f'{zoom}_{tx}_{ty}.png'; url=f'https://tile.openstreetmap.org/{zoom}/{tx}/{ty}.png'; urls.append(url)
            if not file.exists():
                response=session.get(url,timeout=30); response.raise_for_status(); file.write_bytes(response.content)
            mosaic.paste(Image.open(file).convert('RGB'),((tx-tx0)*256,(ty-ty0)*256))
    source=np.asarray(mosaic); dst=np.zeros((1000,1000,3),np.uint8)
    source_transform=from_bounds(tx0*span-world,world-(ty1+1)*span,(tx1+1)*span-world,world-ty0*span,source.shape[1],source.shape[0])
    for band in range(3):
        reproject(source[:,:,band],dst[:,:,band],src_transform=source_transform,src_crs='EPSG:3857',dst_transform=from_bounds(extent[0],extent[2],extent[1],extent[3],1000,1000),dst_crs='EPSG:25832',resampling=Resampling.bilinear)
    s=pd.read_csv(R/'01/survey_legacy_selected.csv').drop_duplicates('ObjektIDObjektID')
    s=s[s.xcoord.between(*extent[:2])&s.ycoord.between(*extent[2:])]
    fig,ax=plt.subplots(figsize=(8,8),layout='constrained'); ax.imshow(dst,extent=extent,origin='upper')
    ax.vlines(np.arange(extent[0],extent[1]+1,10),extent[2],extent[3],color='k',lw=.2,alpha=.16)
    ax.hlines(np.arange(extent[2],extent[3]+1,10),extent[0],extent[1],color='k',lw=.2,alpha=.16)
    ax.scatter(s.xcoord,s.ycoord,s=8,c='#b82538',edgecolors='white',linewidths=.2,label='Boreholes')
    axis(ax); ax.set(xlim=extent[:2],ylim=extent[2:]); ax.legend(loc='upper left')
    ax.text(.99,.01,'© OpenStreetMap contributors',transform=ax.transAxes,ha='right',fontsize=8,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    save(fig,'01_grid_boreholes')
    (O/'data/basemap.json').write_text(json.dumps({'source':'https://www.openstreetmap.org/copyright','tiles':urls,'map_crs':'EPSG:25832','zoom':zoom},indent=2))

if a.part in ['sections','vtk']:
    fields={k:np.load(R/f'04/{k}.npy').reshape(nz,ny,nx) for k in ['s_type','entropy','var']}
    tier={'entropy':classes(fields['entropy'],.3,.7),'var':classes(fields['var'],.2,.4)}
    stats={'entropy_max_below_ground':float(fields['entropy'][mask].max()),'entropy_above_one_below_ground':int((fields['entropy'][mask]>1).sum()),'entropy_classes':'<0.3; 0.3–0.7 inclusive; >0.7 (including values >1, unchanged legacy entropy)','variability_classes':'<0.2; 0.2–0.4 inclusive; >0.4','colors':{'Clay':'#9955bb','Sand':'#f28e2b','Quaternary':'#ffdf00'}}
    (O/'data/uncertainty_classes.json').write_text(json.dumps(stats,indent=2),encoding='utf8')
    if a.part=='sections':
        if a.x_cut is None or a.y_cut is None: p.error('sections require --x-cut and --y-cut in EPSG:25832')
        if not(extent[0]<=a.x_cut<extent[1] and extent[2]<=a.y_cut<extent[3]): p.error('section coordinates outside domain')
        # A requested grid edge uses the cell immediately to its positive side.
        ix=int((a.x_cut-extent[0])/10); iy=int((a.y_cut-extent[2])/10)
        cat=ListedColormap(list(stats['colors'].values())); levels=ListedColormap(['#70b7a4','#f5d76e','#c85a54'])
        cuts=[(np.s_[:,:,ix],y,dem[:,ix],f'X = {a.x_cut:.0f} m','Y [m]'),(np.s_[:,iy,:],x,dem[iy,:],f'Y = {a.y_cut:.0f} m','X [m]')]
        for key,label,bins in [('entropy','Entropy',['<0.3','0.3–0.7','>0.7']),('var','Variability',['<0.2','0.2–0.4','>0.4'])]:
            fig,axs=plt.subplots(2,3,figsize=(16,8),layout='constrained')
            for row,(sl,h,surface,title,hl) in enumerate(cuts):
                for col,(v,name) in enumerate([(fields['s_type'],'Most likely class'),(fields[key],label),(tier[key],label+' class')]):
                    ax=axs[row,col]; opts=dict(cmap=cat if col==0 else levels,norm=BoundaryNorm([.5,1.5,2.5,3.5],3)) if col!=1 else dict(cmap='magma',vmin=0,vmax=float(fields[key][mask].max()) if key=='entropy' else 2/3)
                    im=ax.imshow(np.where(mask[sl],v[sl],np.nan),origin='lower',extent=[h[0]-5,h[-1]+5,z[0]-.25,z[-1]+.25],aspect='auto',interpolation='nearest',**opts)
                    ax.plot(h,surface,c='.2',lw=1,label='Digital elevation model'); ax.set_ylim(z[0]-.25,np.ceil(surface.max()+2)); axis(ax,hl,'Elevation [m]')
                    ax.set_title(name+' · '+title,fontsize=10)
                    if col==0: ax.legend(fontsize=8,loc='lower left')
                    bar=fig.colorbar(im,ax=ax,shrink=.8,pad=.02)
                    if col!=1: bar.set_ticks([1,2,3],labels=list(stats['colors']) if col==0 else bins)
            save(fig,'13_sections_entropy' if key=='entropy' else '14_sections_variability')
        (O/'data/sections.json').write_text(json.dumps({'requested_x':a.x_cut,'requested_y':a.y_cut,'sampled_x_center':float(x[ix]),'sampled_y_center':float(y[iy]),'rotation':0},indent=2))
    else:
        import pyvista as pv
        from vtk import vtkUnstructuredGridWriter
        mesh=pv.read(a.repo/'vtk/model_gok_exponential.vtk'); ids=mesh.cell_data['grid_cell_id']
        for key in ['N300','N300_depth','support_300','confidence_300','priority_300']:
            if key in mesh.cell_data: del mesh.cell_data[key]
        for key in ['entropy','var']: mesh.cell_data[key+'_class']=tier[key].ravel()[ids]
        mesh.cell_data['DEM']=mesh.cell_data.pop('GOK') if 'GOK' in mesh.cell_data else mesh.cell_data['DEM']
        out=O.parent/'model_gok_exponential.vtk'; writer=vtkUnstructuredGridWriter(); writer.SetFileName(str(out)); writer.SetInputData(mesh); writer.SetFileVersion(42); writer.SetFileTypeToBinary(); assert writer.Write()==1
        check=pv.read(out); assert check.n_cells==int(mask.sum())
        for key in fields: np.testing.assert_array_equal(check.cell_data[key],fields[key].ravel()[ids].astype(check.cell_data[key].dtype))
        assert np.all(z[ids//(nx*ny)]+.25<=dem.ravel()[ids%(nx*ny)])
        print('VTK verified:',check.n_cells, 'cells;',out.stat().st_size,'bytes')
