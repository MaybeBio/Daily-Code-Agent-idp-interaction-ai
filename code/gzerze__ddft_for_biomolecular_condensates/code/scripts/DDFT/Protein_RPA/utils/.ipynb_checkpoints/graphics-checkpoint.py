#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modules that are used to define the:

    1. Convert set of images into a GIF + MP4
    2. Generate 3-D images using MAYAVI (no support for cluster currently)
"""

import os
import moviepy.editor as mp
#from mayavi import mlab
#mlab.options.offscreen = True
import numpy as np
import subprocess

def key_func(x):
    """
    Input   =   List of PNG image file paths (x)

    Output  =   List of filenames without extensions (to sort)
    """
    return int(x.split('_')[-1].rstrip('.png'))

def key_func_t(x):
    """
    Input   =   List of TXT image file paths (x)

    Output  =   List of filenames without extensions (to sort)
    """
    return int(x.split('_')[-1].rstrip('.txt'))

def image_list(dir_name,file_names,start_char):
    """
    **Input**
    
        dir_name    =   Path to parent directory of images
        file_names  =   list of sorted image file names
        start_char  =   label to pick (P for protein, R for RNA, c if chi images are plotted)

    Output  =   List of complete path to particular sorted set of images
    """
    images = [];
    for filename in file_names:
        if filename[0]==start_char:
            images.append((dir_name +  filename))

    return(images)

def save_movie(dir_name,out_dir,movie_name ='movie',duration=0.1):
    """
    Function generates 2D movies (gif & MP4) from output images

    **Input parameters**

        -   dir_name    =   path to directory with image files
        -   out_dir    =   path to directory of movie files
        -   movie_name  =   prefix for output movie file (default = movie)
        -   duration    =   time between 2 frames (default = 0.1)
    """
    file_names = list((fn for fn in os.listdir(dir_name) if fn.endswith('.png')))
    file_names = sorted(file_names,key=key_func)

    clip = mp.ImageSequenceClip(image_list(dir_name,file_names,'P'),fps=1/duration);
    clip.write_gif(out_dir +  movie_name + '_P.gif');

    clip = mp.ImageSequenceClip(image_list(dir_name,file_names,'M'),fps=1/duration);
    clip.write_gif(out_dir +  movie_name + '_M.gif');

    clip = mp.ImageSequenceClip(image_list(dir_name,file_names,'D'),fps=1/duration);
    clip.write_gif(out_dir +  movie_name + '_D.gif');

    m_A  =out_dir +'*.gif';
    bash_cmd = 'gifsicle --batch --optimize ' + m_A
    res = subprocess.check_output(['bash','-c',bash_cmd])

    file_names = list((fn for fn in os.listdir(out_dir) if fn.endswith('.gif')))
    for f in file_names:
        clip = mp.VideoFileClip(out_dir +  f)
        clip.write_videofile(out_dir +  f.rstrip('gif') + 'mp4')

def save_movie_1d(dir_name,movie_name ='movie',duration=0.1):
    """
    Function generates 2D movies (gif & MP4) from output images

    **Input parameters**

        -   dir_name    =   path to directory with image files
        -   movie_name  =   prefix for output movie file (default = movie)
        -   duration    =   time between 2 frames (default = 0.1)
    """
    file_names = list((fn for fn in os.listdir(dir_name) if fn.endswith('.png')))
    file_names = sorted(file_names,key=key_func)


    clip = mp.ImageSequenceClip(image_list(dir_name,file_names,'D'),fps=1/duration);
    clip.write_gif(dir_name +  movie_name + '_D.gif');

    clip = mp.ImageSequenceClip(image_list(dir_name,file_names,'M'),fps=1/duration);
    clip.write_gif(dir_name +  movie_name + '_M.gif');

    m_A  =dir_name +'*.gif';
    bash_cmd = 'gifsicle --batch --optimize ' + m_A
    res = subprocess.check_output(['bash','-c',bash_cmd])

    file_names = list((fn for fn in os.listdir(dir_name) if fn.endswith('.gif')))
    for f in file_names:
        clip = mp.VideoFileClip(dir_name +  f)
        clip.write_videofile(dir_name +  f.rstrip('gif') + 'mp4')



def generate_images_3D(dir_name,label_idx=3,N=30,colormap="Blues",vmin=0.0,vmax=0.8,opacity=0.2):
    """
    Function generates 3D images from mesh data

    **No support** for headless nodes (cluster compute nodes for e.g.) and requires *mayavi* to be installed

    **Input parameters**

        -   dir_name    =   path to output_folder
        -   label_idx   =   label of species (3 for protein & 4 for RNA - column number in mesh data)
        -   N           =   length of cubical grid to plot (passed from phase_field function)
        -   colormap    =   colormap to visualize phase-field (default = Blues)
        -   vmin        =   Value of phase-field to align to lower end of color spectrum
        -   vmax        =   Value of phase-field to align to upper end of color spectrum
        -   opacity     =   Transperancy of 3-D rendered image
    """
    mesh_files= (os.listdir(dir_name+'Mesh/'));
    mesh_files = sorted(mesh_files,key=key_func_t)

    data = np.genfromtxt(dir_name+'Mesh/' + mesh_files[0],skip_header=True)
    phi_a = data[:,label_idx].reshape(N,N,N)
#    iso = mlab.contour3d(phi_a,contours=256,opacity=opacity,vmin=vmin,vmax=vmax,colormap=colormap,line_width=0.2);
    iso = mlab.points3d(phi_a,opacity=opacity,vmin=vmin,vmax=vmax,colormap=colormap,line_width=0.2);

    mlab.outline()
    if label_idx==3:
        mlab.savefig(filename=dir_name + 'Images/P_step_0.png')
    elif label_idx==4:
        mlab.savefig(filename=dir_name + 'Images/R_step_0.png')

    for out_file in mesh_files[1:]:

        idx= out_file.rstrip('.txt').split('_')[-1];
        data = np.genfromtxt(dir_name+'Mesh/' + out_file,skip_header=True)

        phi_a = data[:,label_idx].reshape(N,N,N)

        iso.mlab_source.scalars = phi_a
        if label_idx==3:
            mlab.savefig(filename=dir_name + 'Images/P_step_' + str(idx) + '.png')
        elif label_idx==4:
            mlab.savefig(filename=dir_name + 'Images/R_step_' + str(idx) + '.png')
#    mlab.axes()
    mlab.close()
