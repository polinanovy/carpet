import timeit
start = timeit.default_timer()
import pandas as pd
import urllib.request as ur

# assign data of lists. 
Alerts = {'AlertType': ['IceCubeCASCADE',     # --> 0
                        'IceCubeGoldBronze'], # --> 1
                'URL': ['https://gcn.gsfc.nasa.gov/amon_icecube_cascade_events.html',       # --> 0
                        'https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html']}  # --> 1
 
# Create DataFrame 
AlertSrcData = pd.DataFrame(Alerts) 
print('\n\n')
#print(AlertSrcData); print('\n\n')

df = {} # dictionary to store dataframes - generally better than an array

# https://superfastpython.com/multiprocessing-pool-for-loop/
for i in range(0,len(AlertSrcData)):  #range(0,6): 
    start1 = timeit.default_timer()
    print('Parsing: ',AlertSrcData.AlertType[i],'---',AlertSrcData.URL[i])
    if ur.urlopen(AlertSrcData.URL[i]).getcode() == 200:
        print('\tPing - OK')
        print('\tStart reading html')
        tab = pd.read_html(AlertSrcData.URL[i])
        if len(tab) == 1:
            print('\tReading html - OK (time: ',timeit.default_timer()-start1,')')
            df[i] = tab[0]
            #print(df[i])
    else:
        print('\tPing - FALL')

# Formating of Dataframes
print('\n')

# IceCubeCASCADE
df[0] = df[0].drop(columns=[('OBSERVATION', 'SKYMAP')])
df[0] = df[0].drop([(len(df[0])-5),(len(df[0])-4),(len(df[0])-3),(len(df[0])-2),(len(df[0])-1)])
df[0] = df[0].rename(columns={'RunNum_EventNum':'RunEventNum','Time UT':'TimeUT','delta_T':'DeltaT'})
#print((df[0]))
df[0].to_csv('IceCubeCASCADE.csv', sep=';')
print('\n')

# IceCubeGoldBronze
df[1] = df[1].rename(columns={'RunNum_EventNum':'RunEventNum','Time UT':'TimeUT','RA [deg]':'RA','Dec [deg]':'Dec','Error90 [arcmin]':'Error90','Error50 [arcmin]':'Error50','FAR [#/yr]':'FAR'})
#print((df[1]))
df[1].to_csv('IceCubeGoldBronze.csv', sep=';')
print('\n')


print('Total time: ', timeit.default_timer() - start)
