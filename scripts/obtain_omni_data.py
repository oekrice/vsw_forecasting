#This is a testbed for making the forecast scripts actually nice, and doing it all properly and things. HA, that went well!

import time

t0 = time.time()
import os
import sys
import numpy as np
import csv
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
t1 = time.time()

print('Import 1 time', t1-t0)

t2 = time.time()
print('Import 2 time', t2-t1)
#Just a script to download OMNI measurements and save them in a reasonable format.
#Will eventually look into removing CMEs and other snazzy things.
#But this will require some nuance.

#READ THIS
#____________________________________________
#It's a bit of a bodge, but this script needs to be run twice. Do it with redownload_data = True THEN redownload_data = False. It's just easier that way.

redownload_data = False

omni_fname = './data/shared_data/omni.csv'
cme_fname = './data/shared_data/cme_list.csv'


if os.path.exists(omni_fname) and not redownload_data:
    omni_data = []

    with open(omni_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            omni_data.append([datetime.strptime(row[0], "%Y-%m-%d %H:%M:%S"), float(row[1])])

    omni_data = np.array(omni_data)
    all_dtime_omni = omni_data[:,0]
    all_vsw_omni = omni_data[:,1]

    print('OMNI Data loaded from existing file.')

else:
    print('Downloading OMNI data...')
    path = os.getcwd()

    dtime_min = datetime(2005, 1, 1)
    dtime_max = datetime(2025, 9, 1)
    #dtime_max = datetime(2006, 9, 1)

    data_omni = Hin.get_omni(dtime_min, dtime_max)
    all_dtime_omni = data_omni['datetime']
    all_vsw_omni = data_omni['V'].values

#Run through and add information about whether there was a CME at this time. Can then use that information later on.
if os.path.exists(cme_fname):
    cme_data = []
    with open(cme_fname, "r", encoding="utf-8") as f:
        data = csv.reader(f)
        for row in data:
            try:
                cme_data.append([datetime.strptime(row[1], "%Y/%m/%d %H%M"), datetime.strptime(row[2], "%Y/%m/%d %H%M")])
            except:
                pass
            #cme_data.append([row[1],row[2]])

    #This is just a list of the start and end times of each CME.
    cme_data = np.array(cme_data)

else:
    raise Exception('CME List not found')

omni_data = []
cme_counter = 0
cme_time_window_days = 5 #Will remove the CMEs from either side of this window
persistence_cadence = 27.7  #Will flag times which are 27 days after a CME, and as such shouldn't be considered for persistence metric'

#Run through CMEs and add a flag for the times at which they exist. I thought the other way made more sense but I was evidently wrong
cme_flags = np.zeros((len(all_dtime_omni)))
cme_persist_flags = np.zeros((len(all_dtime_omni)))

min_i = 0

if not redownload_data:
    #Logs the times of the CMES
    for ci in range(len(cme_data)):
        cme_start = cme_data[ci][0] - timedelta(days=cme_time_window_days)
        cme_end = cme_data[ci][1] + timedelta(days=cme_time_window_days)
        iscme = 0
        i = min_i
        go = True
        while go and i < len(all_dtime_omni):
            if cme_start < all_dtime_omni[i] and cme_end > all_dtime_omni[i]:
                cme_flags[i] = 1
                if iscme == 0:
                    min_i = i
                iscme = 1
            if all_dtime_omni[i] > cme_end:
                break
            i += 1

    min_i = 0

    #Logs the times of the CMEs, offset by one solar rotation (into the future)
    for ci in range(len(cme_data)):
        cme_start = cme_data[ci][0] - timedelta(days=cme_time_window_days) + timedelta(days=persistence_cadence)
        cme_end = cme_data[ci][1] + timedelta(days=cme_time_window_days) + timedelta(days=persistence_cadence)
        iscme = 0
        i = min_i
        go = True
        while go and i < len(all_dtime_omni):
            if cme_start < all_dtime_omni[i] and cme_end > all_dtime_omni[i]:
                cme_persist_flags[i] = 1
                if iscme == 0:
                    min_i = i
                iscme = 1
            if all_dtime_omni[i] > cme_end:
                break
            i += 1


for i in range(len(all_dtime_omni)):

    omni_data.append([all_dtime_omni[i], all_vsw_omni[i], cme_flags[i], cme_persist_flags[i]])

omni_data = np.array(omni_data)

if True:
    fig1 = plt.figure(figsize=(12,6))

    #Do a sanity check plot
    cme_mask = omni_data[:,2].astype('bool') #There is a CME at this time
    cme_persist_mask = omni_data[:,3].astype('bool') #There is a CME at this time
    #plt.plot(cme_mask)
    # print(cme_mask)
    plt.plot(omni_data[:,0], omni_data[:,1]*cme_mask, c='red')
    plt.plot(omni_data[:,0], omni_data[:,1]*(1.0-np.logical_or(cme_mask, cme_persist_mask)), c='blue')
    plt.plot(omni_data[:,0], omni_data[:,1]*cme_persist_mask, c='green')

    # plt.plot(omni_data[:,0][1-cme_mask][:24*300], omni_data[:,1][1-cme_mask][:24*300], c='red')
    plt.savefig('./plots/cme_times.png')
    plt.show()
#
omni_fname = './data/shared_data/omni.csv'
with open(omni_fname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(omni_data)




