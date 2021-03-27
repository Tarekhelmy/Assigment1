import numpy as np

Cord=[100.000000,0.126000
    ,99.941610,0.134190
    ,99.766580,0.158700
    ,99.475320,0.199380
    ,99.068500,0.255950
    ,98.547090,0.328040
    ,97.912290,0.415190
    ,97.165590,0.516850
    ,96.308730,0.632380
    ,95.343720,0.761080
    ,94.272800,0.902170
    ,93.098490,1.054850
    ,91.823510,1.218230
    ,90.450850,1.391430
    ,88.983720,1.573510
    ,87.425540,1.763530
    ,85.779950,1.960510
    ,84.050790,2.163470
    ,82.242110,2.371420
    ,80.358130,2.583370
    ,78.403240,2.798280
    ,76.382020,3.015150
    ,74.299170,3.232940
    ,72.159580,3.450580
    ,69.968230,3.667000
    ,67.730250,3.881090
    ,65.450850,4.091740
    ,63.135370,4.297780
    ,60.789210,4.498020
    ,58.417860,4.691240
    ,56.026830,4.876190
    ,53.621740,5.051610
    ,51.208190,5.216200
    ,48.791810,5.368660
    ,46.378260,5.507690
    ,43.973170,5.632000
    ,41.582150,5.740330
    ,39.210790,5.831450
    ,36.864630,5.904190
    ,34.549150,5.957470
    ,32.269760,5.990280
    ,30.031770,6.001720
    ,27.840420,5.991020
    ,25.700830,5.957550
    ,23.617990,5.900810
    ,21.596760,5.820480
    ,19.641870,5.716400
    ,17.757890,5.588560
    ,15.949210,5.437150
    ,14.220050,5.262510
    ,12.574460,5.065130
    ,11.016280,4.845670
    ,9.549150,4.604890
    ,8.176490,4.343710
    ,6.901520,4.063100
    ,5.727200,3.764140
    ,4.656280,3.447920
    ,3.691270,3.115590
    ,2.834410,2.768270
    ,2.087710,2.407060
    ,1.452910,2.033000
    ,0.931490,1.647060
    ,0.524680,1.250110
    ,0.233420,0.842890
    ,0.058390,0.426030
    ,0.000000,0.000000]

# function for finding an orthognal 2d vector to vector a
def perpendicular(a):
    b=np.empty_like(a)
    b[0]=-a[1]
    b[1]=a[0]
    return b
# function to normalize the 2d vector
def normalize(a):
    a=np.array(a)
    return a/np.linalg.norm(a)

xcord=Cord[0::2]
xcord=xcord[0::10]
xcord=xcord[::-1]
ycord=Cord[1::2]
ycord=ycord[0::10]
ycord=ycord[::-1]
freev=[50,0]
orthogonalgradientv=[perpendicular(normalize([(xcord[i+1]-xcord[i]),(ycord[i+1]-ycord[i])])) for i in range(len(xcord)-1)]

# to find the normal velocity, the free stream velocity is dot producted with the orthogonal vector of the panel
normv=[np.dot(freev,gradientv[i]) for i in range(len(gradientv))]
xcordcontroldiff=[((-xcord[i]+xcord[i+1])/2000) for i in range(len(xcord)-1)]
ycordcontroldiff=[((-ycord[i]+ycord[i+1])/2000)for i in range(len(ycord)-1)]
# finding the area requirement of the vortex sheet in order to create an equivilant but opposite velocity to vnorm
gammaareareq=[2*-1*np.pi*(np.sqrt((xcordcontroldiff[i]**2) +(ycordcontroldiff[i]**2)))*normv[i] for i in range(len(xcordcontroldiff))]

print(normv)
print(gammaareareq)
gammainit=0
gamma=[]
gammainit=1000
gamma.append(gammainit)
pytha=[]
for i in range(len(xcord)-1):
    # rearanging the area formula of a trapizoid in order to get the change in gamma required and thus
    # adding it to the initial gamma
    gammainit1=gammainit
    pytha.append(np.sqrt(xcordcontroldiff[i]**2 +ycordcontroldiff[i]**2))
    gammainit1+=(gammaareareq[i] -gammainit*pytha[i])*4/(pytha[i])
    gamma.append(gammainit1)


print(gamma)




