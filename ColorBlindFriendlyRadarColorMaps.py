#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 26.02.26                                                 #
# ColorBlindFriendlyRadarColorMaps.py                                         #
#                                                                             #
# Following the suggestion of Sherman et al. 2024 (Kai Mühlbauer as Coauthor; #
# https://doi.org/10.1175/BAMS-D-23-0056.1) here is a Color Map suitable for  #
# everybody INCLUDING people with Color Vision Deficiencies. This is taken    #
# from https://github.com/openradar/cmweather but made for having directly    #
# discrete color maps.                                                        #
#                                                                             #
# Instead of 'ChaseSpectral' the cmap 'HomeyerRainbow' works as well!         #
# --------------------------------------------------------------------------- #

import matplotlib as mpl
import numpy as np
import wradlib as wrl  # for special colormaps: ChaseSpectral

# --------------------------------------------------------------------------- #
# 14 Colors :                                                                 #
# --------------------------------------------------------------------------- #

# Old MIUB Radar Colors:
colors_radar_old = np.array(
    [[0.00, 1.00, 1.00], [0.00, 0.70, 0.93], [0.00, 0.00, 1.00],  # blues
     [0.50, 1.00, 0.00], [0.40, 0.80, 0.00], [0.27, 0.55, 0.00],  # greens
     [1.00, 1.00, 0.00], [0.80, 0.80, 0.00], [1.00, 0.65, 0.00],  # yellows
     [1.00, 0.27, 0.00], [0.80, 0.22, 0.00], [0.55, 0.15, 0.00],  # reds
     [1.00, 0.00, 1.00], [0.58, 0.44, 0.86]])  # pinks
cmap_radar_old = mpl.colors.ListedColormap(colors_radar_old)

# ChaseSpectral:
colors_radar = mpl.colormaps._cmaps['ChaseSpectral'](np.linspace(0, 1, 15))[1:]
cmap_radar = mpl.colors.ListedColormap(colors_radar)

colors_radar_rho = mpl.colormaps._cmaps['ChaseSpectral'](np.linspace(0, 1, 15))[:-1]
cmap_radar_rho = mpl.colors.ListedColormap(colors_radar_rho)

# --------------------------------------------------------------------------- #
# 15 Colors (White + 14):                                                     #
# --------------------------------------------------------------------------- #

# diameters
colors_radar_dm = mpl.colormaps._cmaps['ChaseSpectral_r'](np.linspace(0, 1, 19))[:-5]
colors_radar_dm = mpl.colormaps._cmaps['ChaseSpectral_r'](np.linspace(0, 1, 19))[1:-4]
cmap_radar_dm = mpl.colors.ListedColormap(
    np.append(np.array([[1,1,1,1],]),colors_radar_dm,axis=0))

# I/L Water  content
# colors_radar_cont = (mpl.colormaps._cmaps['ChaseSpectral_r'](np.linspace(0, 1, 19))[1:-4])
colors_radar_cont = (mpl.colormaps._cmaps['ChaseSpectral_r'](np.linspace(0, 1, 19))[1:-4])
colors_radar_cont=np.delete(colors_radar_cont,-8,axis=0) # remove one orange
cmap_radar_cont = mpl.colors.ListedColormap(
    np.append(np.array([[1,1,1,1],]), colors_radar_cont,axis=0))

# --------------------------------------------------------------------------- #
# 13 Colors (White + 12):                                                     #
# --------------------------------------------------------------------------- #

# number concentration
colors_radar_nt = mpl.colormaps._cmaps['ChaseSpectral_r'](np.linspace(0, 1, 15))[2:-1]
cmap_radar_nt = mpl.colors.ListedColormap(
    np.append(np.array([[1,1,1,1],]), colors_radar_nt,axis=0))

# --------------------------------------------------------------------------- #
# levels/norm for pol. moments and retrievals:                                #
# --------------------------------------------------------------------------- #

# Zh (13 ticks)
levels_zh = [-10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
norm_zh = mpl.colors.BoundaryNorm(levels_zh, len(levels_zh) - 1)
# ZDR (13 ticks)
levels_zdr = [-1, -.1, 0, .1, .2, .3, .4, .5, .6, .8, 1, 2, 3]
norm_zdr = mpl.colors.BoundaryNorm(levels_zdr, len(levels_zdr) - 1)
# KDP (13 ticks)
levels_kdp = [-.5, -.1, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 2, 3]
norm_kdp = mpl.colors.BoundaryNorm(levels_kdp, len(levels_kdp) - 1)
# RHOHV (13 ticks)
levels_rhohv = [.7, .8, .85, .9, .92, .94, .95, .96, .97, .98, .99, .995, .998]
norm_rhohv = mpl.colors.BoundaryNorm(levels_rhohv, len(levels_rhohv) - 1)
# Dm (14 ticks)
levels_dm =  np.linspace(0, 6.5, 14)
norm_dm = mpl.colors.BoundaryNorm(levels_dm, len(levels_dm) - 1)
# I/L Water Cont (be careful of my special '0.01'-tick; 13 ticks)
levels_cont = np.array([0, 0.01, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1])
norm_cont = mpl.colors.BoundaryNorm(levels_cont, len(levels_cont) - 1)
# nt (12 ticks)
levels_nt = np.linspace(-3, 2.5, 12)
norm_nt = mpl.colors.BoundaryNorm(levels_nt, len(levels_nt) - 1)


