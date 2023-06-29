# DRF Spectrum Plotter
# Authors: W. Engelke, AB4EJ, University of Alabama
#           Cuong Nguyen, KC3UAX, The University of Scranton

import matplotlib.pyplot as plt
import numpy as np
import digital_rf as drf
import maidenhead, time, datetime, math
from datetime import timezone

currTime = time.time()

plt.style.use('_mpl-gallery-nogrid')
maidenheadGrid = 'EN91'                     # init maidenheadGrid variable

dataDir = "OBS2023-06-26"
metadata_dir = dataDir + '\ch0\metadata'

do = drf.DigitalRFReader(dataDir)
s, e = do.get_bounds('ch0')

t = dataDir[3:13]
print("Plotting spectrum for " + t)
t = t + "T00:00"
requestTime = datetime.datetime.strptime(t, '%Y-%m-%dT%H:%M')
# these are based on unix time * 10 (for 10 samples/sec)
timestamp = requestTime.replace(tzinfo=timezone.utc).timestamp() * 10

s = int(timestamp)
print("time stamp ",s)

freqList = [0]
theLatitude = 0
theLongitude = 0

# ---------------------------------------
# Get metadata
try:
    dmr = drf.DigitalMetadataReader(metadata_dir)
    first_sample, last_sample = dmr.get_bounds()
    print("metadata bounds are %i to %i" % (first_sample, last_sample))

    start_idx = int(np.uint64(first_sample))
    print('computed start_idx = ',start_idx)

    fields = dmr.get_fields()
    print("Available fields are <%s>" % (str(fields)))

    print("first read - just get one column ")
    data_dict = dmr.read(start_idx, start_idx + 2, "center_frequencies")
    for key in data_dict.keys():
        freqList = data_dict[key]
        print("freq = ",freqList[0])

    data_dict = dmr.read(start_idx, start_idx + 2, "lat")
    for key in data_dict.keys():
        theLatitude = data_dict[key]
        print("Latitude: ",theLatitude)
        
    data_dict = dmr.read(start_idx, start_idx + 2, "long")
    for key in data_dict.keys():
        theLongitude = data_dict[key]
        print("Longitude: ",theLongitude)

    maidenheadGrid = maidenhead.to_maiden(theLatitude, theLongitude)
    
except IOError:
    print("IO Error; metadata not found at " + metadata_dir)

# ----------------------------------------
# Read data
gain = 1.5
hr1 = np.arange(1024, dtype='f')
zeros = np.zeros(1024, dtype='f')

print('Read data... this may take a few minutes...')
offset = 0
# Get data from DRF dataset and build the 2D Q array
for i in range(1439): # 1439 gives 1440 bins
    try:
        data = do.read_vector(s + offset, 1024, 'ch0') 
        h = np.fft.fft(data)  
        k = 0
    # convert fft to waterfall using magnitude
        for j in range(511, 2, -1): # invert bins for plotting (511, 0)
            hr1[k] = (math.sqrt(h[j].real ** 2 + h[j].imag ** 2) ) ** (1/gain)
            k += 1
        for j in range(1023, 512, -1): # initially starts at 1023
            hr1[k] = (math.sqrt(h[j].real ** 2 + h[j].imag ** 2)) ** (1/gain)
            k += 1
        if (offset == 0):
            Q = np.array([hr1])
        else:
            Q = np.append(Q, [hr1], axis=0)
    except IOError: # tried to read DRF data but did not find requsted time slice
        if (offset == 0):
            Q = np.array([zeros])
        else:
            Q = np.append(Q, [zeros], axis = 0) # 
        
    # in narrow case, there are 10 samples/sec, so 600 samples = 1 minute
    offset = offset + 600  #  note overlap of the 1024 bins
    if (i % 100 == 0): # progress indicator, marching dots
        print(".",end='')

    Q[i,1022] = 0 # blank out bottom bin

print("File loaded")

print(time.time()-currTime)

# -----------------------------------------------
# Plot

fig, ax = plt.subplots(2, 1, figsize=(14,8)) # plot size, inches x and y

intensity = 1.5 # multiplier to intensify spectrum display
P = intensity * np.absolute(Q.transpose(1, 0)) # re-orient to show horizontal waterfall
Q1 = np.absolute(Q)

D = np.zeros((1024,1440)) # derivative
D1 = np.zeros((1024, 1440))

M = np.zeros((1440))
amplitudes = np.zeros(1440)
for column in range(0, 1439):  # iterate thru minutes
    maxptr = 0
    maxval = 0
    # here we skip topmost 1 or 2 bins as they contiain only the DC offet
    for bin in range(1, 1020): # iterate thru bins of this minute
        if Q1[column, bin] > maxval:
            maxval = Q1[ column, bin ]
            maxptr = bin  # remember which bin had highest amplitude
        D[bin, column] = 5 * (P[bin, column] - P[bin-1, column])
        D1[bin, column] = 5 * (D[bin, column] - D[bin-1, column])

    M[column] = maxptr  # remember max bin for each column
    amplitudes[column] = maxval

avgAmplitude = np.average(amplitudes, weights=None)
stdDevFreq = np.std(M[0:840])
print("First 14 hours STD DEV = ",stdDevFreq)

ax[0].text(40,40,'std.dev. = ' + str(stdDevFreq),color='white')

# original cmap was gist_ncar
ax[0].imshow(P, cmap =plt.get_cmap('viridis'), interpolation ='nearest', vmin=0,vmax=(avgAmplitude*1.5))

plt.subplot(211)

x_points = []
y_points = []

for x in range(1439):
    x_points.append(x)
    y_points.append(M[x])
plt.scatter(x_points, y_points, c='red', marker = "+", s=1) # s is is marker size

ax[0].set_aspect('auto')

plt.yticks(np.arange(0,1100,100), labels=['5','4','3','2','1','0','-1','-2','-3','-4','-5'])
plt.xticks(np.arange(0,1560, 120), labels=['00','02','04','06','08','10','12','14','16','18','20','22','24'])
ax[0].set_ylabel('Hz')
ax[0].set_xlabel('Hours, UTC')
ax[0].set_title('Grape Narrow Spectrum, Freq. = ' + str(freqList[0]) + " MHz, " + t + ' , Lat. '
                + '{:6.2f}'.format(theLatitude) + ", Long. " + '{:6.2f}'.format(theLongitude) + ' (Grid'
                + maidenheadGrid + ')')
plt.subplot(212)
plt.margins(x=0)
plt.autoscale(enable=True,axis='y')
ax[1].plot(amplitudes)

plt.xticks(np.arange(0,1560, 120), labels=['00','02','04','06','08','10','12','14','16','18','20','22','24'])
ax[1].set_ylabel('Amplitude, uncalibrated units')
ax[1].set_xlabel('Hours, UTC')
ax[1].set_title('Peak Amplitude by Minute' )
fig.tight_layout()
print("Average amplitude: ",avgAmplitude)

print(time.time()-currTime)

plt.savefig('images\DRF_plots\\' + t + '_' + maidenheadGrid + '.jpg')
plt.show()