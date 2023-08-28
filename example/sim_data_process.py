#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 10:17:32 2023

@author: ostocker
"""

"""
PyCrown - Fast raster-based individual tree segmentation for LiDAR data
-----------------------------------------------------------------------
Copyright: 2018, Jan Zörner
Licence: GNU GPLv3
"""
import os
import sys
sys.path.append("/home/ostocker/travaux/OSToolBox")
import OSToolBox as ost
from pathlib import Path
from datetime import datetime
from pycrown import PyCrown
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from PIL import Image
import numpy as np
from math import floor
import pandas as pd
import laspy
import plyfile
import matplotlib.pyplot as plt

def ply2las(points, to_extract, output_path):
    data_pts=[]
    to_extract=[to_extract] if not isinstance(to_extract, list) else to_extract
    for i in to_extract:
        idx=np.where(points[:,-1]==i)
        data_pts.append(points[idx])
    tree_points=np.vstack(data_pts)
    try:
        xmin = np.floor(np.min(tree_points[:,0]))
    except:
        return 0
    ymin = np.min(tree_points[:,1])
    zmin = np.min(tree_points[:,2])
    xmax = np.max(tree_points[:,0])
    ymax = np.max(tree_points[:,1])
    zmax = np.max(tree_points[:,2])
    
    x = tree_points[:,0]
    y = tree_points[:,1]
    z = tree_points[:,2]
    intensity = tree_points[:,3]
    classif = tree_points[:,4]
    
    output_path = Path(output_path)
    
    header = laspy.header.LasHeader()
    outfile = laspy.LasData(header)
    outfile.header.minimum = [xmin, ymin, zmin]
    outfile.header.maximum = [xmax, ymax, zmax]
    outfile.x = x
    outfile.y = y
    outfile.z = z
    outfile.intensity = intensity
    outfile.classif = classif
    outfile.write(output_path)
    return outfile
            

def create_dem(point_cloud, offset, maximum, cell_size, out):
    min_coords = [0.0, 0.0, 0.0]
    max_coords = np.array(maximum) - np.array(offset)#[np.max(point_cloud.x), np.max(point_cloud.y)]
    
    # Create a meshgrid for interpolation
    x_grid = np.arange(min_coords[0], max_coords[0], cell_size)
    y_grid = np.arange(min_coords[1], max_coords[1], cell_size)
    
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)

    # Flatten the grid to create a list of points for querying
    query_points = np.vstack((x_grid.flatten(), y_grid.flatten())).T

    # Build a k-d tree from the ground points for efficient nearest-neighbor search
    x = point_cloud.x; y = point_cloud.y; z = point_cloud.z
    np_point_cloud = np.c_[(x, y, z)]

    tree = cKDTree(np_point_cloud[:,:2])
    
    # Set the number of nearest neighbors to use in the IDW interpolation (you can adjust this)
    num_neighbors = 20

    # Calculate the distances and indices of the nearest neighbors for each grid point
    distances, indices = tree.query(query_points, k=num_neighbors)
    
    # Perform IDW interpolation for each grid point
    interpolated_elevations = np.zeros(query_points.shape[0])
    for i in range(query_points.shape[0]):
        dists = distances[i]<0.66
        inds = indices[i][dists]
        if inds.shape[0] != 0:
            interpolated_elevations[i] = np.max(np_point_cloud[inds, 2])
        else:
            interpolated_elevations[i] = np.max(np_point_cloud[indices[i][0], 2])
        
    # Reshape the interpolated elevations back to a 2D grid
    dem = interpolated_elevations.reshape(x_grid.shape)
    #dem = np.flip(dem,0)
    im = Image.fromarray(dem)
    im.save(Path(out))

    return dem
    
    """
    # Interpolate elevation values using scipy's griddata
    PC_xy = np.c_[point_cloud.x, point_cloud.y]
    PC_z = np.array(point_cloud.z)
    dem = griddata(PC_xy, PC_z, (x, y), method='linear')
    
    im = Image.fromarray(dem)
    im.save(Path(out))
    
    return dem
    """


if __name__ == '__main__':
    cell_size = 1
    TSTART = datetime.now()
    increment = -1
    for folder in ost.getDirBySubstr("/media/ostocker/Storage/data/act1/set2/datasets/","plyfiles"):
        for file in ost.getFileByExt(folder+"/train","ply"):
            # try:
            increment += 1
            # Read the PLY file and extract the X, Y, Z coordinates
            ply_data = plyfile.PlyData.read(file)
            x = ply_data['vertex']['x']
            if x.shape[0] < 1000:
                continue
            y = ply_data['vertex']['y']
            z = ply_data['vertex']['z']
            intensity = np.zeros(x.shape)
            classif = ply_data['vertex']['class']
            points = np.c_[x,y,z,intensity,classif]
            
            to_extract = [0,1,2]
            F_LAS_ALL = ply2las(points, to_extract, '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/all.las')
            to_extract = [1,2]
            F_LAS_POINTS = ply2las(points, to_extract, '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/POINTS.las')
            to_extract = 0
            F_LAS_GROUND = ply2las(points, to_extract, '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/ground.las')
            
            if F_LAS_GROUND == 0:
                print("NO GROUND POINTS")
                print("=======================================")
                continue
            
            F_LAS_GROUND.classif += 30-F_LAS_GROUND.classif   # Put classification label of ground equal to -1
    
            out = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/DSM.tif'
            dsm = create_dem(F_LAS_ALL, F_LAS_ALL.header.offset, F_LAS_ALL.header.maximum, cell_size, out)
            out = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/DTM.tif'
            dtm = create_dem(F_LAS_GROUND, F_LAS_ALL.header.offset, F_LAS_ALL.header.maximum, cell_size, out)
            
            if dtm.shape != dsm.shape:
                print("dtm ", dtm.shape)
                print("dsm ", dsm.shape)
                print(file_name)
            
            chm = dsm - dtm
            
            im = Image.fromarray(chm)
            im.save(Path('/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/CHM.tif'))

            F_CHM = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/CHM.tif'
            F_DTM = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/DTM.tif'
            F_DSM = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/DSM.tif'
            F_LAS = '/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output/POINTS.las'

            PC = PyCrown(F_CHM, F_DTM, F_DSM, F_LAS, outpath='/home/ostocker/travaux/unsupervised_boreal/segmentation_tree/output')

            # Cut off edges
            # PC.clip_data_to_bbox((1802200, 1802400, 5467250, 5467450))

            # Smooth CHM with 5m median filter
            PC.filter_chm(2, ws_in_pixels=True)

            # Tree Detection with local maximum filter
            returned = PC.tree_detection(PC.chm, ws=3, ws_in_pixels=True, hmin=2.)
            if returned == 0:
                continue

            # Clip trees to bounding box (no trees on image edge)
            # original extent: 1802140, 1802418, 5467295, 5467490
            # PC.clip_trees_to_bbox(bbox=(1802150, 1802408, 5467305, 5467480))
            # PC.clip_trees_to_bbox(bbox=(1802160, 1802400, 5467315, 5467470))
            PC.clip_trees_to_bbox(inbuf=1)  # inward buffer of 11 metre

            # PARAMETERS DE PY CROWN : A ADAPTER 
            ###############################################
            # Crown Delineation
            PC.crown_delineation(algorithm='dalponteCIRC_numba', th_tree=2., th_seed=0.4, th_crown=0.5, max_crown=10.)

            # Correct tree tops on steep terrain
            PC.correct_tree_tops()

            # Calculate tree height and elevation
            PC.get_tree_height_elevation(loc='top')
            PC.get_tree_height_elevation(loc='top_cor')

            # Screen small trees
            returned = PC.screen_small_trees(hmin=5., loc='top')
            if returned == 0:
                continue
                
            # Convert raster crowns to polygons
            PC.crowns_to_polys_raster()
            name = "trees_increment_{0:0>6}".format(increment)
            last_folder = name[-10:-7]
            PC.crowns_to_polys_smooth(last_folder, name, F_LAS_GROUND, store_ply=True)

            # Check that all geometries are valid
            returned = PC.quality_control()
            if returned == 0:
                continue
            


            # Export results
            PC.export_raster(PC.chm, PC.outpath / 'chm.tif', 'CHM')
            #PC.export_tree_locations(loc='top')
            #PC.export_tree_locations(loc='top_cor')
            #PC.export_tree_crowns(crowntype='crown_poly_raster')
            #PC.export_tree_crowns(crowntype='crown_poly_smooth')

            print(f"Number of trees detected: {len(PC.trees)}")
            print()
            print("=======================================")
            print()
            # except:
            #     print(folder)
            #     print(file)
            #     print(increment)
            #     continue
            break
        break
    TEND = datetime.now()
    print(f'Processing time: {TEND-TSTART} [HH:MM:SS]')