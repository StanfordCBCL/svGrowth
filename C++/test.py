import vessel
import pickle

def saveVessel(vess):
    with open('vessel.pickle', 'wb') as file:
        pickle.dump(vess,file)
    return

def loadVessel():
    with open('vessel.pickle', 'rb') as file:
        vess = pickle.load(file)
    return vess


A = vessel.Vessel(numLen=4,length=1)
A.initializeVessel()
#A.runFluidIteration()

# with open('vessel.pickle', 'wb') as file:
#     pickle.dump(A,file)

#with open('vessel.pickle', 'rb') as file:
#    A = pickle.load(file)

A.runFluidSolidIteration()
#A.runFluidSolidIteration()

saveVessel(A)

import pdb; pdb.set_trace()
