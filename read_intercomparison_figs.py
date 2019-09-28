import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys

"""
Code to read .jpeg files containing the plots from the BOMEX intercomparison
paper by Siebesma et al 2003, and convert them to timeseries with numbers. Then
output the numbers to a text file.

Example syntax:

    "python read_intercomparison_figs.py <filename> xy <output_name>"

Note:
    Do not put file extensions, this code assumes .jpeg input and .txt output
"""

# click to form a box and remove the figure labels
def select_box(black, grey, my_title = 'Remove Label'):
    """
    Plots the cropped panel, then user clicks around the figure label to return 
    the coordinates of the corners of a polygon surrounding the label. Values
    within that polygon are reset to 0.0 and the panel is returned.
    """
    corners = []
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    ax.contourf(black, levels = [0.5, 1.1], colors = ['k'])
    ax.contourf(grey, levels = [0.5, 1.1], colors = ['grey'], alpha = 0.33)
    ax.set_title(my_title)
    plt.pause(0.1)
    corners.append(plt.ginput(4))
    
    # convert the corners to ints in x and y
    x_corners = [int(corner[0]) for corner in corners[0]]
    y_corners = [int(corner[1]) for corner in corners[0]]
    ax.plot([x_corners[idx] for idx in range(-1, len(x_corners))], [y_corners[idx] for idx in range(-1, len(y_corners))], 'r--')
    plt.pause(1)
    plt.close('all')
    return x_corners, y_corners

# The grid dots in the original figure make this retrieval `spiky' attempt to remove those spikes
def remove_spikes(X):
    """
    First removes spikes based on local gradients, then uses a median filter on 
    the data with a 25-point window to remove remaining spikes
    """
    # remove spikes based on local gradients, if its changing unrealistically
    for idx in range(1, len(X)-1):
        if abs(X[idx] - X[idx-1]) > 3*abs(X[idx+1] - X[idx-1]):
            # in-fill the spikes with the mean of the neighbors
            X[idx] = 0.5*(X[idx+1] + X[idx-1])
    X_temp = X*1.
    for idx in range(12, len(X) - 12):
        X[idx] = np.median(X_temp[(idx-12):(idx+12)])
    
    for idx in range(4, len(X) - 4):
        X[idx] = np.mean(X[(idx-4):(idx+4)])
    
    plt.close('all')
    
    return X

def read_jpeg_fig(fig_name, orientation, out_name):
    """
    Uses the Pillow module (PIL) to import a figure as an array.
    1. Assumes that there is onely one panel to the figure
    2. Only greyscale line plots with three colors (i.e. white, grey, and black)
       are supported.
        i. There is support for reading a range around the main line.
    """
    
    ### Read the figure as an array
    fig_in = Image.open(fig_name + '.jpeg', 'r')
    j_size, i_size = fig_in.size
    fig_pix = fig_in.load() 
    # PIL reads as a value from 0 - 255 where 0 = black, and 255 = white
    fig_arr = np.array([[fig_pix[j,i] for j in range(j_size)] for i in range(i_size)])
    
    # Expect the figure to be in greyscale with three colors 
    # (i.e. white, grey, and black)
    
    # Jpeg format has a blur around the colours, use a histogram to get a range
    # of pixel values that correspond to each color
    hist = plt.hist(fig_arr.flatten(), bins = 255)
    plt.close('all')
    
    # Define the maxima from the histogram:
    maxima = [hist[1][idx] for idx in range(len(hist[0])) if ([hist[0][idx-1] < hist[0][idx] if (idx - 1 >= 0) else True][0]) and ([hist[0][idx+1] < hist[0][idx] if (idx + 1 < len(hist[0])) else True][0]) and (hist[0][idx] > 0.001*j_size*i_size)]
    
    # Force black to be zero, white to be 255
    maxima[0] = 0
    maxima[-1] = 255
    
    print maxima
    # Assume a range 25 either side of the maxima and those are the values of the
    # white, grey, and black pixels
    tolerance = 42.5
    # Want the location of the black lines, i.e. the lowest value from the maxima
    # Since the maxima are in ascending order we can do:
    black_arr = np.flipud(np.where((maxima[0] - 2*tolerance <= fig_arr)*(fig_arr <= maxima[0] + 2*tolerance), 1.0, 0.0))
    grey_arr  = np.flipud(np.where((maxima[1] - tolerance <= fig_arr)*(fig_arr <= maxima[1] + tolerance), 1.0, 0.0))
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    ax.contourf(black_arr, levels = [-0.5, 0.5, 1.5], colors = ['None', 'k'])
    ax.contourf(grey_arr, levels = [-0.5, 0.5, 1.5], colors = ['None', 'grey'])
    ax.set_title('Full Image')
    plt.show()
    plt.close('all')
    
    # Plot the black and grey areas
    # Use this plot to manually crop to the figure panel
    l_crop = True
    while l_crop:
        x_crop, y_crop = select_box(black_arr, grey_arr, my_title = 'Select Panel')
        
        # x_min, x_max, y_min, y_max
        x_min, x_max, y_min, y_max = [np.min(x_crop), np.max(x_crop), np.min(y_crop), np.max(y_crop)]
        black_arr = black_arr[y_min:y_max, x_min:x_max]
        grey_arr  = grey_arr[y_min:y_max, x_min:x_max]
        l_crop    = bool(int(raw_input('Do you need to crop again? (No = 0, Yes = 1)')))
    
    # crop to the correct area, assume the curve for the mean is in black, and
    # and assume there is a shaded area in grey for the spread about the mean
    panel_mean   = black_arr*1.
    panel_spread = grey_arr*1.
    
    # If there are figure lettering (e.g. 'a)', 'b)'), remove them
    l_remove_label = bool(int(raw_input('Is there a figure label to remove? (0 = No, 1 = Yes)')))
    while l_remove_label:
        label_corners_x, label_corners_y = select_box(panel_mean, panel_spread)
        panel_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
        panel_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
        l_remove_label = bool(int(raw_input('Is there anything else to remove? (0 = No, 1 = Yes)')))
    
    x_range_n = float(raw_input('What is the lower limit on the x-range?'))
    x_range_x = float(raw_input('What is the upper limit on the x-range?'))
    y_range_n = float(raw_input('What is the lower limit on the y-range?'))
    y_range_x = float(raw_input('What is the upper limit on the y-range?'))
    
    # Get the x-coordinates of the pixels
    x = np.linspace(x_range_n, x_range_x, x_max - x_min)
    
    # Get the y-coordinates of the pixels for each panel
    y = np.linspace(y_range_n, y_range_x, y_max - y_min)
    
    if orientation == 'xy':
        # Interpolate from pixel index to the numbers using the above x- and y-coordinates
        Y_mean  = np.array([y[np.where(panel_mean[:,i] == 1)[0]].mean() if 1 in panel_mean[:,i] else 0.0 for i in range(panel_mean.shape[1])])
        Y_lower = np.array([y[np.where(panel_spread[:,i] == 1)[0]].min() if 1 in panel_spread[:,i] else 0.0 for i in range(panel_spread.shape[1])])
        Y_upper = np.array([y[np.where(panel_spread[:,i] == 1)[0]].max() if 1 in panel_spread[:,i] else 0.0 for i in range(panel_spread.shape[1])])
        
        # If the data is spiky (e.g. due to noise or gridlines), smooth it out
        Y_mean  = remove_spikes(Y_mean)
        Y_lower = remove_spikes(Y_lower)
        Y_upper = remove_spikes(Y_upper)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, Y_mean, 'k', lw = 2)
        ax.fill_between(x, Y_lower, Y_upper, facecolor = 'grey', edgecolor = 'None', alpha = 0.33)
        ax.set_title('Cropped Image With processed data')
        plt.show()
        plt.close('all')
        
        # Save to a txt file
        with open(out_name + '.txt', 'w') as txt_file:
            txt_file.write(out_name + '\n')
            txt_file.write('X : \n')
            [txt_file.write(str(x[i]) + ',') if (i < (len(x) - 1)) else txt_file.write(str(x[i])) for i in range(len(x))]
            txt_file.write('\n')
            txt_file.write('Lower Bounds : \n')
            [txt_file.write(str(Y_lower[i]) + ',') if (i < (len(Y_lower) - 1)) else txt_file.write(str(Y_lower[i])) for i in range(len(Y_lower))]
            txt_file.write('\n')
            txt_file.write('Mean : \n')
            [txt_file.write(str(Y_mean[i]) + ',') if (i < (len(Y_mean) - 1)) else txt_file.write(str(Y_mean[i])) for i in range(len(Y_mean))]
            txt_file.write('\n')
            txt_file.write('Upper Bounds : \n')
            [txt_file.write(str(Y_upper[i]) + ',') if (i < (len(Y_upper) - 1)) else txt_file.write(str(Y_upper[i])) for i in range(len(Y_upper))]
            txt_file.write('\n')
        
    elif orientation == 'yx':
        # Interpolate from pixel index to the numbers using the above x- and y-coordinates
        X_mean  = np.array([x[np.where(panel_mean[j,:] == 1)[0]].mean() if 1 in panel_mean[j,:] else 0.0 for j in range(panel_mean.shape[0])])
        X_lower = np.array([x[np.where(panel_spread[j,:] == 1)[0]].min() if 1 in panel_spread[j,:] else 0.0 for j in range(panel_spread.shape[0])])
        X_upper = np.array([x[np.where(panel_spread[j,:] == 1)[0]].max() if 1 in panel_spread[j,:] else 0.0 for j in range(panel_spread.shape[0])])
        
        # If the data is spiky (e.g. due to noise or gridlines), smooth it out
        X_mean  = remove_spikes(X_mean)
        X_lower = remove_spikes(X_lower)
        X_upper = remove_spikes(X_upper)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(X_mean, y, 'k', lw = 2)
        ax.fill_betweenx(y, X_lower, X_upper, facecolor = 'grey', edgecolor = 'None', alpha = 0.33)
        ax.set_title('Cropped Image With processed data')
        plt.show()
        plt.close('all')
        
        # Save to a txt file
        with open(out_name + '.txt', 'w') as txt_file:
            txt_file.write(out_name + '\n')
            txt_file.write('Y : \n')
            [txt_file.write(str(y[i]) + ',') if (i < (len(y) - 1)) else txt_file.write(str(y[i])) for i in range(len(y))]
            txt_file.write('\n')
            txt_file.write('Lower Bounds : \n')
            [txt_file.write(str(X_lower[i]) + ',') if (i < (len(X_lower) - 1)) else txt_file.write(str(X_lower[i])) for i in range(len(X_lower))]
            txt_file.write('\n')
            txt_file.write('Mean : \n')
            [txt_file.write(str(X_mean[i]) + ',') if (i < (len(X_mean) - 1)) else txt_file.write(str(X_mean[i])) for i in range(len(X_mean))]
            txt_file.write('\n')
            txt_file.write('Upper Bounds : \n')
            [txt_file.write(str(X_upper[i]) + ',') if (i < (len(X_upper) - 1)) else txt_file.write(str(X_upper[i])) for i in range(len(X_upper))]
            txt_file.write('\n')

read_jpeg_fig(sys.argv[1], sys.argv[2], sys.argv[3])

"""
###### Figure 2 ######
fig2 = Image.open('Siebesma03-f02.jpeg', 'r')
j_size, i_size = fig2.size
fig2_pix = fig2.load()
fig2_arr = np.array([[fig2_pix[j,i] for j in range(j_size)] for i in range(i_size)])

# This figure has three panels
# These figures all have three colours (white, grey, and black)
# Jpeg format has a blur around the colours, use a histogram to get a range of 
# pixel values which correspond to each colour
hist = plt.hist(fig2_arr.flatten(), bins = 255)
plt.show()

# find local maxima
maxima = [hist[1][idx] for idx in range(len(hist[0])) if ([hist[0][idx-1] < hist[0][idx] if (idx - 1 >= 0) else True][0]) and ([hist[0][idx+1] < hist[0][idx] if (idx + 1 < len(hist[0])) else True][0]) and (hist[0][idx] > 0.001*j_size*i_size)]

# Assume a range 20 either side of the maxima and those are the values of the
# white, grey, and black pixels

# Want the location of the black lines, i.e. the lowest value from the maxima
# Since the maxima are in ascending order we can do:
arrs = [np.where((max_val - 20 <= fig2_arr)*(fig2_arr <= max_val + 20), 1.0, 0.0) for max_val in maxima]
black_arr, grey_arr, white_arr = [np.flipud(arr) for arr in arrs]

# Create a plot, use it manually define coordinate ranges corresponding to the
# figure panel
plt.contourf(black_arr, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
plt.title('Find the corners of the panels')
plt.show()

# x_min, x_max, y_min, y_max
x_min, x_max, y_min_a, y_max_a = [85, 630, 583, 798] # all the panels have the same x range
y_min_b, y_max_b = [321, 536]
y_min_c, y_max_c = [59, 274]

# test
panel_A_mean   = black_arr[y_min_a:y_max_a, x_min:x_max]
panel_A_spread = grey_arr[y_min_a:y_max_a, x_min:x_max]
panel_B_mean   = black_arr[y_min_b:y_max_b, x_min:x_max]
panel_B_spread = grey_arr[y_min_b:y_max_b, x_min:x_max]
panel_C_mean   = black_arr[y_min_c:y_max_c, x_min:x_max]
panel_C_spread = grey_arr[y_min_c:y_max_c, x_min:x_max]

fig = plt.figure()
axa = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
axa.contourf(panel_A_mean, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
axb = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
axb.contourf(panel_B_mean, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
axc = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
axc.contourf(panel_C_mean, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
axa.set_title('Test all')
plt.show()


label_corners_x, label_corners_y = select_box(panel_A_mean)
panel_A_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_A_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_B_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_B_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_C_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_C_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0

x_range = [0.0, 360.0]   # all panels have the same x range
# Get the x-coordinates of the pixels
x = np.linspace(x_range[0], x_range[1], x_max - x_min)

# Get the y-coordinates of the pixels for each panel
y_range_a = [0, 0.2] # read off the figure
y_a = np.linspace(y_range_a[0], y_range_a[1], y_max_a - y_min_a)

y_range_b = [0, 20]
y_b = np.linspace(y_range_b[0], y_range_b[1], y_max_b - y_min_b)

y_range_c = [0, 750]
y_c = np.linspace(y_range_c[0], y_range_c[1], y_max_c - y_min_c)

# Get the numbers from black lines (the mean of the intercomparison)
# Get the numbers from the grey area (the spread of the intercomparison, +/- 2sigma)
# Panel (a)
Y_a_mean  = np.array([y_a[np.where(panel_A_mean[:,i] == 1)[0]].mean() if 1 in panel_A_mean[:,i] else 0.0 for i in range(panel_A_mean.shape[1])])
Y_a_lower = np.array([y_a[np.where(panel_A_spread[:,i] == 1)[0]].min() if 1 in panel_A_spread[:,i] else 0.0 for i in range(panel_A_spread.shape[1])])
Y_a_upper = np.array([y_a[np.where(panel_A_spread[:,i] == 1)[0]].max() if 1 in panel_A_spread[:,i] else 0.0 for i in range(panel_A_spread.shape[1])])

# Panel (b)
Y_b_mean  = np.array([y_b[np.where(panel_B_mean[:,i] == 1)[0]].mean() if 1 in panel_B_mean[:,i] else 0.0 for i in range(panel_B_mean.shape[1])])
Y_b_lower = np.array([y_b[np.where(panel_B_spread[:,i] == 1)[0]].min() if 1 in panel_B_spread[:,i] else 0.0 for i in range(panel_B_spread.shape[1])])
Y_b_upper = np.array([y_b[np.where(panel_B_spread[:,i] == 1)[0]].max() if 1 in panel_B_spread[:,i] else 0.0 for i in range(panel_B_spread.shape[1])])

# Panel (c)
Y_c_mean  = np.array([y_c[np.where(panel_C_mean[:,i] == 1)[0]].mean() if 1 in panel_C_mean[:,i] else 0.0 for i in range(panel_C_mean.shape[1])])
Y_c_lower = np.array([y_c[np.where(panel_C_spread[:,i] == 1)[0]].min() if 1 in panel_C_spread[:,i] else 0.0 for i in range(panel_C_spread.shape[1])])
Y_c_upper = np.array([y_c[np.where(panel_C_spread[:,i] == 1)[0]].max() if 1 in panel_C_spread[:,i] else 0.0 for i in range(panel_C_spread.shape[1])])

# Save to a txt file
with open('Figure_2_data.txt', 'w') as txt_file:
    txt_file.write('Figure 2 from Siebesma et al 2003\n')
    txt_file.write('Panel A\n')
    txt_file.write('Times : \n')
    [txt_file.write(str(x[i]) + ',') if (i < (len(x) - 1)) else txt_file.write(str(x[i])) for i in range(len(x))]
    txt_file.write('\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(Y_a_lower[i]) + ',') if (i < (len(Y_a_lower) - 1)) else txt_file.write(str(Y_a_lower[i])) for i in range(len(Y_a_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(Y_a_mean[i]) + ',') if (i < (len(Y_a_mean) - 1)) else txt_file.write(str(Y_a_mean[i])) for i in range(len(Y_a_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(Y_a_upper[i]) + ',') if (i < (len(Y_a_upper) - 1)) else txt_file.write(str(Y_a_upper[i])) for i in range(len(Y_a_upper))]
    txt_file.write('\n')
    txt_file.write('Panel B\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(Y_b_lower[i]) + ',') if (i < (len(Y_b_lower) - 1)) else txt_file.write(str(Y_b_lower[i])) for i in range(len(Y_b_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(Y_b_mean[i]) + ',') if (i < (len(Y_b_mean) - 1)) else txt_file.write(str(Y_b_mean[i])) for i in range(len(Y_b_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(Y_b_upper[i]) + ',') if (i < (len(Y_b_upper) - 1)) else txt_file.write(str(Y_b_upper[i])) for i in range(len(Y_b_upper))]
    txt_file.write('\n')
    txt_file.write('Panel C\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(Y_c_lower[i]) + ',') if (i < (len(Y_c_lower) - 1)) else txt_file.write(str(Y_c_lower[i])) for i in range(len(Y_c_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(Y_c_mean[i]) + ',') if (i < (len(Y_c_mean) - 1)) else txt_file.write(str(Y_c_mean[i])) for i in range(len(Y_c_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(Y_c_upper[i]) + ',') if (i < (len(Y_c_upper) - 1)) else txt_file.write(str(Y_c_upper[i])) for i in range(len(Y_c_upper))]
    txt_file.write('\n')

###### Figure 3 ######
# We only really care about panel D here...
fig3 = Image.open('Siebesma03-f03.jpeg', 'r')
j_size, i_size = fig3.size
fig3_pix = fig3.load()
fig3_arr = np.array([[fig3_pix[j,i] for j in range(j_size)] for i in range(i_size)])

# This figure has three panels
# These figures all have three colours (white, grey, and black)
# Jpeg format has a blur around the colours, use a histogram to get a range of 
# pixel values which correspond to each colour
hist = plt.hist(fig3_arr.flatten(), bins = 255)
plt.show()

# find local maxima
maxima = [hist[1][idx] for idx in range(len(hist[0])) if ([hist[0][idx-1] < hist[0][idx] if (idx - 1 >= 0) else True][0]) and ([hist[0][idx+1] < hist[0][idx] if (idx + 1 < len(hist[0])) else True][0]) and (hist[0][idx] > 0.001*j_size*i_size)]
maxima[0] = 0
# Assume a range 85 either side of the maxima and those are the values of the
# white, grey, and black pixels

# Want the location of the black lines, i.e. the lowest value from the maxima
# Since the maxima are in ascending order we can do:
black_arr = np.flipud(np.where((maxima[0] - 85 <= fig3_arr)*(fig3_arr <= maxima[0] + 85), 1.0, 0.0))
grey_arr  = np.flipud(np.where((maxima[1] - 30 <= fig3_arr)*(fig3_arr <= maxima[1] + 30), 1.0, 0.0))

# Create a plot, use it manually define coordinate ranges corresponding to the
# figure panel
plt.contourf(grey_arr, levels = [0.0, 0.5, 1.1], colors = ['w', 'grey'])
plt.title('Find the corners of the panels')
plt.show()

# x_min, x_max, y_min, y_max
x_min, x_max, y_min, y_max = [394, 637, 33, 278] # only doing one panel

# test
panel_mean   = black_arr[y_min:y_max, x_min:x_max]
panel_spread = grey_arr[y_min:y_max, x_min:x_max]

fig = plt.figure()
axa = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
axa.contourf(panel_mean, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
axa.contourf(panel_spread, levels = [0.0, 0.5, 1.1], colors = ['None', 'grey'])
axa.set_title('Test all')
plt.show()

# click to form a box and remove the figure labels
label_corners_x, label_corners_y = select_box(panel_mean)
panel_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0

x_range = [0.0, 0.01]   # all panels have the same x range
# Get the x-coordinates of the pixels
x = np.linspace(x_range[0], x_range[1], x_max - x_min)

# Get the y-coordinates of the pixels for each panel
y_range = [0, 2500] # read off the figure
y = np.linspace(y_range[0], y_range[1], y_max - y_min)

# Get the numbers from black lines (the mean of the intercomparison)
# Get the numbers from the grey area (the spread of the intercomparison, +/- 2sigma)
X_mean  = np.array([x[np.where(panel_mean[j,:] == 1)[0]].mean() if 1 in panel_mean[j,:] else 0.0 for j in range(panel_mean.shape[0])])
X_lower = np.array([x[np.where(panel_spread[j,:] == 1)[0]].min() if 1 in panel_spread[j,:] else 0.0 for j in range(panel_spread.shape[0])])
X_upper = np.array([x[np.where(panel_spread[j,:] == 1)[0]].max() if 1 in panel_spread[j,:] else 0.0 for j in range(panel_spread.shape[0])])

X_mean  = remove_spikes(X_mean)
X_lower = remove_spikes(X_lower)
X_upper = remove_spikes(X_upper)

# Save to a txt file
with open('Figure_3_data.txt', 'w') as txt_file:
    txt_file.write('Figure 3 from Siebesma et al 2003\n')
    txt_file.write('Panel D\n')
    txt_file.write('Height : \n')
    [txt_file.write(str(y[i]) + ',') if (i < (len(y) - 1)) else txt_file.write(str(y[i])) for i in range(len(y))]
    txt_file.write('\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(X_lower[i]) + ',') if (i < (len(X_lower) - 1)) else txt_file.write(str(X_lower[i])) for i in range(len(X_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(X_mean[i]) + ',') if (i < (len(X_mean) - 1)) else txt_file.write(str(X_mean[i])) for i in range(len(X_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(X_upper[i]) + ',') if (i < (len(X_upper) - 1)) else txt_file.write(str(X_upper[i])) for i in range(len(X_upper))]
    txt_file.write('\n')

###### Figure 4 ######
fig4 = Image.open('Siebesma03-f04.jpeg', 'r')
j_size, i_size = fig4.size
fig4_pix = fig4.load()
fig4_arr = np.array([[fig4_pix[j,i] for j in range(j_size)] for i in range(i_size)])

# This figure has three panels
# These figures all have three colours (white, grey, and black)
# Jpeg format has a blur around the colours, use a histogram to get a range of 
# Use the same maxima from figure 3

# Assume a range 85 either side of the maxima and those are the values of the
# white, grey, and black pixels

# Want the location of the black lines, i.e. the lowest value from the maxima
# Since the maxima are in ascending order we can do:
black_arr = np.flipud(np.where((maxima[0] - 85 <= fig4_arr)*(fig4_arr <= maxima[0] + 85), 1.0, 0.0))
grey_arr  = np.flipud(np.where((maxima[1] - 30 <= fig4_arr)*(fig4_arr <= maxima[1] + 30), 1.0, 0.0))

# Create a plot, use it manually define coordinate ranges corresponding to the
# figure panel
plt.contourf(black_arr, levels = [0.0, 0.5, 1.1], colors = ['w', 'k'])
plt.title('Find the corners of the panels')
plt.show()

# x_min, x_max, y_min, y_max
x_min_a, x_max_a, y_min_a, y_max_a = [56, 288, 607, 841]
x_min_b, x_max_b, y_min_b, y_max_b = [350, 580, 607, 841]
#x_min_c, x_max_c, y_min_c, y_max_c = [56, 285, 322, 555]
#x_min_d, x_max_d, y_min_d, y_max_d = [350, 576, 326, 552]
#x_min_e, x_max_e, y_min_e, y_max_e = [56, 285, 38, 270]

# test
panel_A_mean   = black_arr[y_min_a:y_max_a, x_min_a:x_max_a]
panel_A_spread = grey_arr[y_min_a:y_max_a, x_min_a:x_max_a]
panel_B_mean   = black_arr[y_min_b:y_max_b, x_min_b:x_max_b]
panel_B_spread = grey_arr[y_min_b:y_max_b, x_min_b:x_max_b]

# click to form a box and remove the figure labels
label_corners_x, label_corners_y = select_box(panel_A_mean)
panel_A_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_A_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0

label_corners_x, label_corners_y = select_box(panel_B_mean)
panel_B_mean[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0
panel_B_spread[min(label_corners_y):max(label_corners_y), min(label_corners_x):max(label_corners_x)] = 0.0

x_range_a = [0.0, 175.0]
x_range_b = [-40.0, 10.0]
# Get the x-coordinates of the pixels
x_a = np.linspace(x_range_a[0], x_range_a[1], x_max_a - x_min_a)
x_b = np.linspace(x_range_b[0], x_range_b[1], x_max_b - x_min_b)

# Get the y-coordinates of the pixels for each panel
y_range = [0, 2500] # all panels have the same y-range
y_a = np.linspace(y_range[0], y_range[1], y_max_a - y_min_a)
y_b = np.linspace(y_range[0], y_range[1], y_max_b - y_min_b)

# Get the numbers from black lines (the mean of the intercomparison)
# Get the numbers from the grey area (the spread of the intercomparison, +/- 2sigma)
X_A_mean  = np.array([x_a[np.where(panel_A_mean[j,:] == 1)[0]].mean() if 1 in panel_A_mean[j,:] else 0.0 for j in range(panel_A_mean.shape[0])])
X_A_lower = np.array([x_a[np.where(panel_A_spread[j,:] == 1)[0]].min() if 1 in panel_A_spread[j,:] else 0.0 for j in range(panel_A_spread.shape[0])])
X_A_upper = np.array([x_a[np.where(panel_A_spread[j,:] == 1)[0]].max() if 1 in panel_A_spread[j,:] else 0.0 for j in range(panel_A_spread.shape[0])])

X_B_mean  = np.array([x_b[np.where(panel_B_mean[j,:] == 1)[0]].mean() if 1 in panel_B_mean[j,:] else 0.0 for j in range(panel_B_mean.shape[0])])
X_B_lower = np.array([x_b[np.where(panel_B_spread[j,:] == 1)[0]].min() if 1 in panel_B_spread[j,:] else 0.0 for j in range(panel_B_spread.shape[0])])
X_B_upper = np.array([x_b[np.where(panel_B_spread[j,:] == 1)[0]].max() if 1 in panel_B_spread[j,:] else 0.0 for j in range(panel_B_spread.shape[0])])

X_A_mean  = remove_spikes(X_A_mean)
X_A_lower = remove_spikes(X_A_lower)
X_A_upper = remove_spikes(X_A_upper)

X_B_mean  = remove_spikes(X_B_mean)
X_B_lower = remove_spikes(X_B_lower)
X_B_upper = remove_spikes(X_B_upper)

plt.plot(X_A_mean, y_a)
plt.fill_betweenx(y_a, X_A_lower, X_A_upper, facecolor = 'grey', edgecolor = 'None', alpha = 0.33)
plt.show()

plt.plot(X_B_mean, y_b)
plt.fill_betweenx(y_b, X_B_lower, X_B_upper, facecolor = 'grey', edgecolor = 'None', alpha = 0.33)
plt.show()

# Save to a txt file
with open('Figure_4_data.txt', 'w') as txt_file:
    txt_file.write('Figure 4 from Siebesma et al 2003\n')
    txt_file.write('Panel A\n')
    txt_file.write('Height : \n')
    [txt_file.write(str(y_a[i]) + ',') if (i < (len(y_a) - 1)) else txt_file.write(str(y_a[i])) for i in range(len(y_a))]
    txt_file.write('\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(X_A_lower[i]) + ',') if (i < (len(X_A_lower) - 1)) else txt_file.write(str(X_A_lower[i])) for i in range(len(X_A_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(X_A_mean[i]) + ',') if (i < (len(X_A_mean) - 1)) else txt_file.write(str(X_A_mean[i])) for i in range(len(X_A_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(X_A_upper[i]) + ',') if (i < (len(X_A_upper) - 1)) else txt_file.write(str(X_A_upper[i])) for i in range(len(X_A_upper))]
    txt_file.write('\n')
    txt_file.write('Panel B\n')
    txt_file.write('Height : \n')
    [txt_file.write(str(y_b[i]) + ',') if (i < (len(y_b) - 1)) else txt_file.write(str(y_b[i])) for i in range(len(y_b))]
    txt_file.write('\n')
    txt_file.write('Lower Bounds : \n')
    [txt_file.write(str(X_B_lower[i]) + ',') if (i < (len(X_B_lower) - 1)) else txt_file.write(str(X_B_lower[i])) for i in range(len(X_B_lower))]
    txt_file.write('\n')
    txt_file.write('Mean : \n')
    [txt_file.write(str(X_B_mean[i]) + ',') if (i < (len(X_B_mean) - 1)) else txt_file.write(str(X_B_mean[i])) for i in range(len(X_B_mean))]
    txt_file.write('\n')
    txt_file.write('Upper Bounds : \n')
    [txt_file.write(str(X_B_upper[i]) + ',') if (i < (len(X_B_upper) - 1)) else txt_file.write(str(X_B_upper[i])) for i in range(len(X_B_upper))]
    txt_file.write('\n')
"""



