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

import numpy as np
import matplotlib as mpl
# import wradlib as wrl  # for special colormaps: ChaseSpectral

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

# ChaseSpectral ZH ZDR KDP
# colors_radar = mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 15))[1:]
colors_radar = np.array([
    [0.1483253 , 0.12702331, 0.216649 , 1.],
    [0.27216775, 0.23075019, 0.42925388, 1.],
    [0.40089137, 0.34444079, 0.65969808, 1.],
    [0.24987215, 0.56946133, 0.73133195, 1.],
    [0.45168071, 0.74199467, 0.58734849, 1.],
    [0.77096273, 0.87287801, 0.55814971, 1.],
    [0.99129025, 0.98303277, 0.71612057, 1.],
    [0.93607141, 0.76128782, 0.44140871, 1.],
    [0.92727713, 0.49384066, 0.27536666, 1.],
    [0.8102349 , 0.23012616, 0.27361338, 1.],
    [0.63270419, 0.06681815, 0.35554858, 1.],
    [0.87950861, 0.41977631, 0.66841314, 1.],
    [0.62152896, 0.4073525 , 0.68406589, 1.],
    [0.25355039, 0.        , 0.29884244, 1.]])
cmap_radar = mpl.colors.ListedColormap(colors_radar)

# ChaseSpectral CC
# colors_radar_rho = mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 15))[:-1]
colors_radar_rho = np.array(
    [[0.00249159, 0.        , 0.01238498, 1.        ],
     [0.1483253 , 0.12702331, 0.216649  , 1.        ],
     [0.27216775, 0.23075019, 0.42925388, 1.        ],
     [0.40089137, 0.34444079, 0.65969808, 1.        ],
     [0.24987215, 0.56946133, 0.73133195, 1.        ],
     [0.45168071, 0.74199467, 0.58734849, 1.        ],
     [0.77096273, 0.87287801, 0.55814971, 1.        ],
     [0.99129025, 0.98303277, 0.71612057, 1.        ],
     [0.93607141, 0.76128782, 0.44140871, 1.        ],
     [0.92727713, 0.49384066, 0.27536666, 1.        ],
     [0.8102349 , 0.23012616, 0.27361338, 1.        ],
     [0.63270419, 0.06681815, 0.35554858, 1.        ],
     [0.87950861, 0.41977631, 0.66841314, 1.        ],
     [0.62152896, 0.4073525 , 0.68406589, 1.        ]])
cmap_radar_rho = mpl.colors.ListedColormap(colors_radar_rho)

# --------------------------------------------------------------------------- #
# 15 Colors (White + 14):                                                     #
# --------------------------------------------------------------------------- #

# diameters
# colors_radar_dm = mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 19))[ 4:-1]
# colors_radar_dm=np.append(
#     np.array([[1, 1, 1, 1], ]), colors_radar_dm, axis=0)
colors_radar_dm = np.array(
    [[1.        , 1.        , 1.        , 1.        ],
     [0.40708423, 0.36167776, 0.68118068, 1.        ],
     [0.25395644, 0.54662847, 0.74352072, 1.        ],
     [0.34871029, 0.69366826, 0.60511045, 1.        ],
     [0.59993411, 0.80013964, 0.58139943, 1.        ],
     [0.83332996, 0.902003  , 0.56306363, 1.        ],
     [0.99129025, 0.98303277, 0.71612057, 1.        ],
     [0.94044792, 0.8144617 , 0.49306084, 1.        ],
     [0.93122106, 0.61949377, 0.33151531, 1.        ],
     [0.9108918 , 0.39691512, 0.24722056, 1.        ],
     [0.79016854, 0.20211936, 0.27671434, 1.        ],
     [0.61108033, 0.01473938, 0.30621619, 1.        ],
     [0.84366522, 0.26347059, 0.58019299, 1.        ],
     [0.86625306, 0.60778072, 0.80957898, 1.        ],
     [0.55909005, 0.30067979, 0.61448309, 1.        ]])
cmap_radar_dm = mpl.colors.ListedColormap(colors_radar_dm)

# I/L Water  content
# colors_radar_cont = (mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 19))[4:-1])
# colors_radar_cont = np.delete(colors_radar_cont, -8,
#                               axis=0)  # remove one orange
# colors_radar_cont=np.append(
#     np.array([[1, 1, 1, 1], ]), colors_radar_cont, axis=0)
colors_radar_cont = np.array(
    [[1.        , 1.        , 1.        , 1.        ],
     [0.40708423, 0.36167776, 0.68118068, 1.        ],
     [0.25395644, 0.54662847, 0.74352072, 1.        ],
     [0.34871029, 0.69366826, 0.60511045, 1.        ],
     [0.59993411, 0.80013964, 0.58139943, 1.        ],
     [0.83332996, 0.902003  , 0.56306363, 1.        ],
     [0.99129025, 0.98303277, 0.71612057, 1.        ],
     [0.93122106, 0.61949377, 0.33151531, 1.        ],
     [0.9108918 , 0.39691512, 0.24722056, 1.        ],
     [0.79016854, 0.20211936, 0.27671434, 1.        ],
     [0.61108033, 0.01473938, 0.30621619, 1.        ],
     [0.84366522, 0.26347059, 0.58019299, 1.        ],
     [0.86625306, 0.60778072, 0.80957898, 1.        ],
     [0.55909005, 0.30067979, 0.61448309, 1.        ]])
cmap_radar_cont = mpl.colors.ListedColormap(colors_radar_cont)

# --------------------------------------------------------------------------- #
# 13 Colors (White + 12):                                                     #
# --------------------------------------------------------------------------- #

# number concentration
# colors_radar_nt = mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 15))[3:]
# colors_radar_nt=np.append(
#     np.array([[1, 1, 1, 1], ]), colors_radar_nt, axis=0)
colors_radar_nt = np.array(
    [[1.        , 1.        , 1.        , 1.        ],
     [0.40089137, 0.34444079, 0.65969808, 1.        ],
     [0.24987215, 0.56946133, 0.73133195, 1.        ],
     [0.45168071, 0.74199467, 0.58734849, 1.        ],
     [0.77096273, 0.87287801, 0.55814971, 1.        ],
     [0.99129025, 0.98303277, 0.71612057, 1.        ],
     [0.93607141, 0.76128782, 0.44140871, 1.        ],
     [0.92727713, 0.49384066, 0.27536666, 1.        ],
     [0.8102349 , 0.23012616, 0.27361338, 1.        ],
     [0.63270419, 0.06681815, 0.35554858, 1.        ],
     [0.87950861, 0.41977631, 0.66841314, 1.        ],
     [0.62152896, 0.4073525 , 0.68406589, 1.        ],
     [0.25355039, 0.        , 0.29884244, 1.        ]])
cmap_radar_nt = mpl.colors.ListedColormap(colors_radar_nt)

# dm huge
# colors_radar_dm2 = mpl.colormaps._cmaps['ChaseSpectral_r'](
#     np.linspace(0, 1, 19))[2:-4]
# colors_radar_dm2=np.append(
#     np.array([[1, 1, 1, 1], ]), colors_radar_dm2, axis=0)
colors_radar_dm2 = np.array(
    [[1.        , 1.        , 1.        , 1.        ],
     [0.86625306, 0.60778072, 0.80957898, 1.        ],
     [0.84366522, 0.26347059, 0.58019299, 1.        ],
     [0.61108033, 0.01473938, 0.30621619, 1.        ],
     [0.79016854, 0.20211936, 0.27671434, 1.        ],
     [0.9108918 , 0.39691512, 0.24722056, 1.        ],
     [0.93122106, 0.61949377, 0.33151531, 1.        ],
     [0.94044792, 0.8144617 , 0.49306084, 1.        ],
     [0.99685864, 0.99483133, 0.72638202, 1.        ],
     [0.83332996, 0.902003  , 0.56306363, 1.        ],
     [0.59993411, 0.80013964, 0.58139943, 1.        ],
     [0.34871029, 0.69366826, 0.60511045, 1.        ],
     [0.25395644, 0.54662847, 0.74352072, 1.        ],
     [0.40708423, 0.36167776, 0.68118068, 1.        ]])
cmap_radar_dm2 = mpl.colors.ListedColormap(colors_radar_dm2)

# --------------------------------------------------------------------------- #
# 9 Colors (White + 8):                                                     #
# --------------------------------------------------------------------------- #

# number concentration
# colors_radar_nt2 = mpl.colormaps._cmaps['ChaseSpectral'](
#     np.linspace(0, 1, 11))[1:-1]
# colors_radar_nt2=np.append(
#     np.array([[1, 1, 1, 1], ]), colors_radar_nt2, axis=0)
colors_radar_nt2 = np.array(
    [[1.        , 1.        , 1.        , 1.        ],
     [0.19542234, 0.1662855 , 0.29644516, 1.        ],
     [0.3814893 , 0.32358123, 0.62220743, 1.        ],
     [0.26333873, 0.60250045, 0.70176995, 1.        ],
     [0.65276865, 0.8219142 , 0.57547435, 1.        ],
     [0.99129025, 0.98303277, 0.71612057, 1.        ],
     [0.93240992, 0.66345522, 0.36064926, 1.        ],
     [0.8386172 , 0.27196469, 0.26637978, 1.        ],
     [0.68748034, 0.11546206, 0.42190222, 1.        ],
     [0.7672592 , 0.57823069, 0.79554752, 1.        ]])
cmap_radar_nt2 = mpl.colors.ListedColormap(colors_radar_nt2)

# --------------------------------------------------------------------------- #
# levels/norm for pol. moments and retrievals:                                #
# --------------------------------------------------------------------------- #

# ZH (13 ticks)
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
levels_dm = np.linspace(0, 6.5, 14)
norm_dm = mpl.colors.BoundaryNorm(levels_dm, len(levels_dm) - 1)
# Dm huge(13 ticks)
levels_dm2 = np.linspace(0, 18, 13)
norm_dm2 = mpl.colors.BoundaryNorm(levels_dm2, len(levels_dm2) - 1)
# Dm rain (14 ticks)
levels_dm_r = np.linspace(0, 2.6, 14)
norm_dm_r = mpl.colors.BoundaryNorm(levels_dm_r, len(levels_dm_r) - 1)

# I/L Water Cont (be careful of my special '0.01'-tick; 13 ticks)
levels_cont = np.array(
    [0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1])
norm_cont = mpl.colors.BoundaryNorm(levels_cont, len(levels_cont) - 1)

# nt (12 ticks)
levels_nt = np.linspace(-3, 2.5, 12)
norm_nt = mpl.colors.BoundaryNorm(levels_nt, len(levels_nt) - 1)
# nt (9 ticks)
levels_nt2 = np.linspace(-7, 1., 9)
norm_nt2 = mpl.colors.BoundaryNorm(levels_nt2, len(levels_nt2) - 1)
