"""Deterministic fictional boreholes for a redistributable workflow demonstration.

Coordinates overlap Munich numerically; values are generated, NOT observations.
This is a software fixture, never engineering evidence.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.affinity import rotate


def generate(destination):
    destination = Path(destination)/"data"/"bohrungen_epsg25832_shp"
    if destination.exists(): raise FileExistsError(destination)
    destination.mkdir(parents=True)
    rng = np.random.RandomState(20260912)
    geometry, collars, layers = [], [], []
    for j in range(10):
        for i in range(10):
            bh = f"SYNTHETIC_{j:02d}_{i:02d}"
            x, y = 688555.+100*i, 5335505.+100*j
            elevation = 550.+.001*(x-689000)-.001*(y-5336000)
            base_depth = 20.+2*np.sin(i/2.)+1.2*np.cos(j/2.)+rng.normal(0,.5)
            geometry.append(rotate(Point(x,y),-25.67,origin=[688000,5332000]))
            collars.append(dict(ObjektID=bh,Ansatzhoeh=elevation,Endteufe=100.,Gemeinde="SYNTHETIC",Bohrjahr=2026))
            boundaries = [0.,base_depth,*np.arange(30.,101.,10.)]
            for k,(top,bottom) in enumerate(zip(boundaries[:-1],boundaries[1:])):
                din = "G" if k == 0 else ("S" if (i+j+k)%3 else "T")
                layers.append(dict(ObjektID=bh,Obergrenz=top,Untergrenz=bottom,DIN=din,
                                   PetBez={"G":"synthetic gravel","S":"synthetic sand","T":"synthetic clay"}[din],Strati="synthetic"))
    pd.DataFrame(collars).to_csv(destination/"bohrungen_stammdaten.csv",sep=";",decimal=",",index=False)
    pd.DataFrame(layers).to_csv(destination/"bohrungen_schichten.csv",sep=";",decimal=",",index=False)
    gpd.GeoDataFrame({"ObjektID":[r["ObjektID"] for r in collars]},geometry=geometry,crs="EPSG:25832").to_file(destination/"bohrungen.shp",encoding="UTF-8")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    generate(p.parse_args().output)
