"""
Obtained from http://wiki.scipy.org/Cookbook/Matplotlib/ColormapTransformations

"""
import numpy as np
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
import colorsys

"""
Following is a custom colormap for radar plots. It is a rainbow colormap 
but without much green, and saturation and luminance are forced to vary
montonically and evenly from bottom to top.
"""
rgb = []                     
rgb.append([0.000, 0.167, 1.000])
rgb.append([0.100, 0.400, 1.000])
rgb.append([0.200, 0.600, 1.000])
rgb.append([0.400, 0.800, 1.000])
rgb.append([0.600, 0.933, 1.000])
rgb.append([0.800, 1.000, 1.000])
rgb.append([1.000, 1.000, 0.800]) 
rgb.append([1.000, 0.933, 0.600]) 
rgb.append([1.000, 0.800, 0.400]) 
rgb.append([1.000, 0.600, 0.200])  
rgb.append([1.000, 0.400, 0.100]) 
rgb.append([1.000, 0.167, 0.000])  
rgbarray = np.asarray(rgb) 
Lcoeff = 0.65 - np.arange(12) * 0.025
for i, row in enumerate(rgbarray):
    h, l, s = colorsys.rgb_to_hls(row[0], row[1], row[2])
    l *= Lcoeff[i] / l 
    s = 0.70 + i * 0.025
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    #print 'RGB', r, g, b, 'HLS', h, l ,s
    rgbarray[i] = np.asarray([r, g, b])
nrows = rgbarray.shape[0]
xvals = np.linspace(0,1,nrows)
#Now create dictionary of list of tuples (x,y0,y1), x=R, G, or B
cdict = {}
cdict['red'] = [tuple([xvals[x],rgbarray[x,0],rgbarray[x,0]]) for x in np.arange(nrows)]
cdict['green'] = [tuple([xvals[x],rgbarray[x,1],rgbarray[x,1]]) for x in np.arange(nrows)]
cdict['blue'] = [tuple([xvals[x],rgbarray[x,2],rgbarray[x,2]]) for x in np.arange(nrows)]
Lang12 = LinearSegmentedColormap('testcmap', cdict)

def cmap_map(function,cmap):

    """ 
    Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     
    points in a colormap.

    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         
        step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(map( reduced_cmap, step_list))
    new_LUT = np.array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def lighten_cmap(cmap):
    return cmap_map(lambda x: x/2+0.5, cmap)


