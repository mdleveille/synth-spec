#make lists of g and T values so we can read in the models.

glist = ["600","625","650","675","700","725","750","775","800","825",
         "850","875","900"]

garray100 = np.array((glist), dtype=np.float)
garray = [g/100.0 for g in garray100]
Tlist = ["11000","12000","13000","14000","16000","18000","20000",
         "22000","24000","26000","28000","30000","35000"]
Tarray = np.array((Tlist), dtype=np.float)

modinfo = np.array((glist, Tlist))
modflux1 = np.zeros((len(glist), len(Tlist), len(dataf1c)))
modflux2 = np.zeros((len(glist), len(Tlist), len(dataf2c)))

#read in the models for BS1
for i in xrange(0, len(glist)):
    for j in xrange(0, len(Tlist)):
        file = "../Natalie/da" + str(Tlist[j]) +"_" + \
                       str(glist[i]) +".dk.rebin1.lsf.bin30"
        modw1, modf1, moderr1 = np.loadtxt(file, unpack=True)
        dummy = ma.masked_array(modf1, mask=dataw1.mask) 
        modflux1[i][j] = 1.0e-8*dummy.compressed()
                #the model fluxes are per cm. the above
                #line changes this to per A.