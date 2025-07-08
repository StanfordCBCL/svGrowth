from utils import *

# This class provides functionality to running an FSG simulation
class Vessel():

    def __init__(self, **kwargs):

        self.vesselType = kwargs.get("vesselType", "cylinder")
        self.radius = kwargs.get("radius", 1.0)
        self.thickness = kwargs.get("thickness", 0.1)
        self.length = kwargs.get("length", 10.0)
        self.vesselSegmentations = kwargs.get("vesselSegmentations", "")
        self.vesselReference = None
        self.vesselSolid = None
        self.vesselFluid = None
        self.solidResult = None
        self.fluidResult = None

        self.nG = 8
        self.max_days = kwargs.get("gnrMaxDays", 720)
        self.prefix = kwargs.get("prefix", "./")
        self.numRad = kwargs.get("numRad", 4)
        self.numCirc = kwargs.get("numCirc", 12)
        self.numLen = kwargs.get("numLen", 12)
        self.aneurysm = kwargs.get("aneurysm", 0)
        self.tevg = kwargs.get("tevg", 0)
        self.timeStep = 0
        self.timeIter = 0
        self.omega = 0.1
        self.residual = 1.0
        self.simulationInputDirectory = kwargs.get("simulationInputDirectory", "FolderSimulationInputFiles")
        self.simulationExecutable = kwargs.get("simulationExecutable","~/svFSI-build/svFSI-build/mysvfsi")
        self.vesselName = kwargs.get("vesselName", "vessel")
        self.simulationType = kwargs.get("simulationType", "mm")
        self.resultNum = kwargs.get("resultNum",1000)
        self.resultDir = kwargs.get("resultDir","results")
        self.fluidDir = kwargs.get("fluidDir","fluid-results")
        self.outputDir = kwargs.get("outputDir","Outputs")
        self.estimateFluids = kwargs.get("estimateFluids", True)
        self.predictMethod = kwargs.get("predictMethod", "aitken")
        self.gnr_step_size = kwargs.get("gnrStepSize", 1.0)
    

    def incrementTimestep():
        self.timeStep+=1
        return

    def incrementIteration():
        self.timeIter+=1
        return


    def runFluidIteration(self):
        self.updateFluid()
        self.saveFluid()
        self.runFluid()
        self.updateFluidResults()

        return


    def runFluidSolidIteration(self):
        self.updateSolid()
        self.saveSolid()
        self.updateFluid()
        self.appendFluidResult()
        self.saveFluid()
        self.runFluidSolid()
        self.updateFluidSolidResults()
        self.appendSolidResult()
        self.appendIterfaceResult()
        self.updateMaterial()
        self.updateReference()

        self.timeIter+=1
        if checkTol:
            self.timeStep+=1

        return


    def runSolidIteration(self):
        self.updateSolid()
        self.saveSolid()
        self.runSolid()
        self.updateSolidResults()
        self.updateMaterial()
        self.updateReference()

        self.timeIter+=1
        if checkTol:
            self.timeStep+=1
            
        return

    def checkTol(self):
        return


    def runSolid(self):
        os.system("mpiexec " + self.simulationExecutable + " " + self.simulationInputDirectory + "/solid_"+self.simulationType+".mfs")
        return

    def runFluid(self):
        os.system("mpiexec " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_fluid.mfs")
        return

    def runFluidSolid(self):
        os.system("mpiexec " + self.simulationExecutable + " " + self.simulationInputDirectory + "/input_"+self.simulationType+".mfs")
        return

    def initializeVessel(self):
        if self.vesselType == "cylinder":
            self.vesselReference = self.initializeCylinder()
            self.updateSolid()
        else:
            print("TODO: Implement segmentation initializer")
        return

    def updateMaterial(self):


        numCells = self.vesselReference.GetNumberOfCells()
        input_array = []

        for q in range(numCells):
            sigma_inv = self.vesselReference.GetCellData().GetArray('inv_curr').GetTuple1(q)
            wss = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            simulate = int(self.vesselReference.GetCellData().GetArray('Simulate').GetTuple1(q))
            aneurysmValue = self.vesselReference.GetCellData().GetArray('aneurysmValue').GetTuple1(q)
            tevgValue = self.vesselReference.GetCellData().GetArray('tevgValue').GetTuple1(q)

            #Rotate into GnR membrane frame
            e_r = self.vesselReference.GetCellData().GetArray('e_r').GetTuple(q)
            e_t = self.vesselReference.GetCellData().GetArray('e_t').GetTuple(q)
            e_z = self.vesselReference.GetCellData().GetArray('e_z').GetTuple(q)
            Q = np.array((e_r,e_t,e_z))


            defGrad_mem = []
            #Gauss point values
            for p in range(self.nG):

                defGrad = self.vesselReference.GetCellData().GetArray('defGrad').GetTuple(q)
                defGrad_g = defGrad[p*9:(p+1)*9]
                defGrad_g = np.reshape(defGrad_g, (3,3))
                defGrad_mem_g = np.matmul(Q,np.matmul(defGrad_g, np.transpose(Q)))
                defGrad_mem_g = np.reshape(defGrad_mem_g,9)

                #Queue process
                input_array.append([q, p,self.outputDir, self.timeStep,self.timeIter,simulate,self.max_days,self.gnr_step_size,sigma_inv,wss,aneurysmValue,tevgValue,\
                                   defGrad_mem_g[0],defGrad_mem_g[1],defGrad_mem_g[2],defGrad_mem_g[3],\
                                   defGrad_mem_g[4],defGrad_mem_g[5],defGrad_mem_g[6],defGrad_mem_g[7],defGrad_mem_g[8]])

                defGrad_mem = np.hstack((defGrad_mem,defGrad_mem_g))

            self.vesselReference.GetCellData().GetArray('defGrad_mem').SetTuple(q,defGrad_mem)

        input_file = open("input_array.dat","wb")
        pickle.dump(input_array, input_file)
        input_file.close()

        print("Running points...")
        os.system("mpiexec python utils_run_vessel.py")

        print("Parsing points...")
        pool = Pool()
        # Begin processing
        output_array = []
        for output in tqdm(pool.imap(parsePoint,input_array), total=len(input_array)):
            output_array.append(output)
        pool.close()
        pool.join()
        pool.terminate()

        for q in range(numCells):
            stiffness_mem_g = []
            sigma_mem_g = []
            J_di_g = []
            for p in range(self.nG):
                J_di, stiffness, sigma, wss, sigma_inv = output_array[q*self.nG + p]

                stiffness_mem_g = stiffness_mem_g + stiffness.tolist()
                sigma_mem_g = sigma_mem_g + sigma.tolist()
                J_di_g = J_di_g + [J_di]

            self.vesselReference.GetCellData().GetArray('stiffness_mem').SetTuple(q,stiffness_mem_g)
            self.vesselReference.GetCellData().GetArray('sigma_mem').SetTuple(q,sigma_mem_g)
            self.vesselReference.GetCellData().GetArray('J_di').SetTuple(q,J_di_g)

        return

    def updateSolid(self):
        self.vesselSolid = self.vesselReference.warp_by_vector("displacements")
        arrayNames = self.vesselReference.array_names
        for name in arrayNames:
            if name not in ["GlobalNodeID", "varWallProps", "GlobalElementID", "InnerRegionID", "OuterRegionID", "DistalRegionID", "ProximalRegionID", "StructureID", "Pressure"]:
                if name in self.vesselSolid.point_data:
                    self.vesselSolid.point_data.remove(name)
                if name in self.vesselSolid.cell_data:
                    self.vesselSolid.cell_data.remove(name)
        return

    def updateFluid(self):
        surf = self.vesselSolid.extract_surface()
        inner = thresholdModel(surf, 'InnerRegionID',0.5,1.5)
        self.vesselFluid = self.generateFluidMesh(inner)

    def saveSolid(self):
        vol = self.vesselSolid

        os.makedirs(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces', exist_ok=True)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-complete.mesh.vtu',vol)
        surf = vol.extract_surface()
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-complete.mesh.vtp',surf)
        outer = thresholdModel(surf, 'OuterRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_outer.vtp',outer)
        distal = thresholdModel(surf, 'DistalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_inlet.vtp',distal)
        proximal = thresholdModel(surf, 'ProximalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_outlet.vtp',proximal)
        inner = thresholdModel(surf, 'InnerRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/solid-mesh-complete/mesh-surfaces/wall_inner.vtp',inner)

        return

    def saveFluid(self):
        fluidVol = self.vesselFluid

        os.makedirs(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces', exist_ok=True)
        # Add global node id to inner volume
        numPtsVol = fluidVol.GetNumberOfPoints()
        numCellsVol = fluidVol.GetNumberOfCells()

        # Need to order ids to match globalNodeID order
        fluidVol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPtsVol,numPtsVol).astype(int),name="GlobalNodeID"))
        fluidVol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCellsVol,numCellsVol).astype(int),name="GlobalElementID"))

        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-complete.mesh.vtu',fluidVol)
        innerVolSurf = pv.wrap(fluidVol).extract_surface()
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-complete.mesh.vtp',innerVolSurf)

        innerProx = thresholdModel(innerVolSurf, 'ProximalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_inlet.vtp',pv.wrap(innerProx))
        innerDist = thresholdModel(innerVolSurf, 'DistalRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_outlet.vtp',innerDist)
        innerWall = thresholdModel(innerVolSurf, 'OuterRegionID',0.5,1.5)
        save_data(self.prefix + 'mesh/lumen-mesh-complete/mesh-surfaces/lumen_wall.vtp',innerWall)

        return

    def saveReference(self, saveIter=False):
        if saveIter:
            save_data(self.prefix + 'meshIterations/mesh_' + str(self.timeStep) + '_' + str(self.timeIter) + '.vtu', exist_ok=True)
        else:
            save_data(self.prefix + 'meshResults/mesh_' + str(self.timeStep) + '.vtu', exist_ok=True)
        return


    def updateFluidSolidResults(self):
        resultdata = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        self.solidResult = thresholdModel(resultdata,'Domain_ID', 1.5, 2.5, extract=False, cell=True)
        self.fluidResult = thresholdModel(resultdata,'Domain_ID', 0.5, 1.5, extract=False, cell=True)

        return

    def updateSolidResults(self):
        self.solidResult = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        return

    def updateFluidResults(self):
        self.fluidResult = read_data(self.resultDir+'/result_'+str(self.resultNum)+'.vtu', file_format="vtu")
        return


    def appendIterfaceResult(self):
        """
        TODO
        """
        numCells = self.vesselReference.GetNumberOfCells()

        fluidSurface = self.fluidResult.extract_surface()
        pointLocatorFluid = vtk.vtkPointLocator()
        pointLocatorFluid.SetDataSet(fluidSurface)
        pointLocatorFluid.BuildLocator()

        for q in range(numCells):
            fluidStressId = int(self.vesselReference.GetCellData().GetArray('fluidStressQueryID').GetTuple1(q))
            cell = self.vesselReference.GetCell(fluidStressId)
            cellPts = cell.GetPointIds()

            # For each cell, use trapezoidal integration to compute WSS
            cellWSS = 0.0
            numberOfPoints = 0
            for p in range(cellPts.GetNumberOfIds()):
                ptId = cellPts.GetId(p)
                if self.vesselReference.GetPointData().GetArray('InnerRegionID').GetTuple1(ptId) == 1:
                    fluidId = int(pointLocatorFluid.FindClosestPoint(self.vesselReference.GetPoint(ptId)))
                    pointWSS = fluidSurface.GetPointData().GetArray('WSS').GetTuple3(fluidId)
                    pointPressure = fluidSurface.GetPointData().GetArray('Pressure').GetTuple1(fluidId)
                    self.vesselReference.GetPointData().GetArray('Pressure').SetTuple1(ptId,pointPressure)
                    cellWSS += np.linalg.norm(pointWSS)
                    numberOfPoints += 1
            #Average WSS in cell
            cellWSS *= 1/float(numberOfPoints)

            wss_prev = self.vesselReference.GetCellData().GetArray('wss_curr').GetTuple1(q)
            self.vesselReference.GetCellData().GetArray('wss_prev').SetTuple1(q,wss_prev)
            self.vesselReference.GetCellData().GetArray('wss_curr').SetTuple1(q,cellWSS)

        return


    def estimateIterfaceResult(self):
        """
        TODO
        """
        numPts = self.vesselReference.GetNumberOfPoints()

        for q in range(numPts):
            fluidStressId = int(volumedata.GetCellData().GetArray('fluidStressQueryID').GetTuple1(q))
            cell = fluidSurface.GetCell(fluidStressId)
            cellPts = cell.GetPointIds()

            # For each cell, use trapezoidal integration to compute WSS
            cellWSS = 0.0
            numberOfPoints = 0
            for p in range(cellPts.GetNumberOfIds()):
                ptId = cellPts.GetId(p)
                if self.vesselReference.GetPointData().GetArray('InnerRegionID').GetTuple1(ptId) == 1:
                    fluidCoordinate =  self.vesselReference.GetPoint(ptId)
                    pointWSS = 4*0.04*20.0/(np.pi*(np.linalg.norm(fluidCoordinate[0:2])**3.0))
                    cellWSS += pointWSS
                    numberOfPoints += 1
            #Average WSS in cell
            cellWSS *= 1/float(numberOfPoints)

            wss_prev = volumedata.GetCellData().GetArray('wss_curr').GetTuple1(q)
            volumedata.GetCellData().GetArray('wss_prev').SetTuple1(q,wss_prev)
            volumedata.GetCellData().GetArray('wss_curr').SetTuple1(q,cellWSS)

        return



    def appendSolidResult(self):
        pointLocatorSolid = vtk.vtkPointLocator()
        pointLocatorSolid.SetDataSet(self.solidResult)
        pointLocatorSolid.BuildLocator()

        numPts = self.vesselReference.GetNumberOfPoints()
        numCells = self.vesselReference.GetNumberOfCells()

        for q in range(numPts):
            originalCoordinate = np.array(self.vesselReference.GetPoint(q))
            displacement_prev = np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(q))
            currentCoordinate = originalCoordinate + displacement_prev

            pointIdSolid = pointLocatorSolid.FindClosestPoint(currentCoordinate)
            solidCoordinate = np.array(self.solidResult.GetPoint(pointIdSolid))
            solidDispacement = np.array(self.solidResult.GetPointData().GetArray("Displacement").GetTuple3(pointIdSolid))

            displacement_curr = solidCoordinate + solidDispacement - originalCoordinate

            residual_curr = displacement_curr - displacement_prev
            self.vesselReference.GetPointData().GetArray("residual_curr").SetTuple(q, residual_curr)

        if self.predictMethod == "none":
            self.omega = 1.0
        elif self.predictMethod == "aitken":
            rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr")).flatten()
            rprev = np.array(self.vesselReference.GetPointData().GetArray("residual_prev")).flatten()
            diff = rcurr - rprev
            self.omega = -self.omega*np.dot(rprev,diff)/np.dot(diff,diff)

        # Calculate cauchy green tensor
        for q in range(numPts):
            rcurr = np.array(self.vesselReference.GetPointData().GetArray("residual_curr").GetTuple3(q))
            dcurr = np.array(self.vesselReference.GetPointData().GetArray("displacements").GetTuple3(q))
            displacement = dcurr + self.omega *rcurr

            self.vesselReference.GetPointData().GetArray("residual_prev").SetTuple(q, rcurr)
            self.vesselReference.GetPointData().GetArray("displacements").SetTuple(q, displacement)

            self.vesselReference = computeGaussValues(self.vesselReference,"displacements")


        # Get stress invariant
        for q in range(numCells):
            cell = self.solidResult.GetCell(q)
            cellPts = cell.GetPointIds()
            # For each cell, use trapezoidal integration to compute sigma
            cellSigma = 0.0
            for r in range(cellPts.GetNumberOfIds()):
                ptId = cellPts.GetId(r)
                pointSigma = self.solidResult.GetPointData().GetArray('Cauchy_stress').GetTuple6(ptId)           
                cellSigma += (pointSigma[0]+pointSigma[1]+pointSigma[2])
            cellSigma *= 1/float(cellPts.GetNumberOfIds())

            inv_prev = self.vesselReference.GetCellData().GetArray("inv_curr").GetTuple1(q)
            self.vesselReference.GetCellData().GetArray("inv_prev").SetTuple1(q, inv_prev)
            self.vesselReference.GetCellData().GetArray("inv_curr").SetTuple1(q, cellSigma)

        return

    def appendFluidResult(self):

        numPoints = self.vesselFluid.GetNumberOfPoints()

        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet(self.fluidResult)
        pointLocator.BuildLocator()

        velocity_array = np.zeros((numPoints,3))
        pressure_array = np.zeros((numPoints,1))

        for q in range(numPoints):
            coordinate = self.vesselFluid.GetPoint(q)
            pointIdSource = pointLocator.FindClosestPoint(coordinate)

            velocity_array[q,:] = self.fluidResult.GetPointData().GetArray('Velocity').GetTuple(pointIdSource)
            pressure_array[q,:] = self.fluidResult.GetPointData().GetArray('Pressure').GetTuple1(pointIdSource)

        self.vesselFluid.GetPointData().AddArray(pv.convert_array(np.array(velocity_array),name="Velocity"))
        self.vesselFluid.GetPointData().AddArray(pv.convert_array(np.array(pressure_array),name="Pressure"))

        return


    def updateReference(self):

        numCells = self.vesselReference.GetNumberOfCells()
        numPoints = self.vesselReference.GetNumberOfPoints()

        array_names = self.vesselReference.array_names

        if "temp_array" in array_names:
            self.vesselReference.rename_array("varWallProps","material_global")
            self.vesselReference.rename_array("temp_array","varWallProps")

        for q in range(numCells):
            e_r = self.vesselReference.GetCellData().GetArray('e_r').GetTuple(q)
            e_t = self.vesselReference.GetCellData().GetArray('e_t').GetTuple(q)
            e_z = self.vesselReference.GetCellData().GetArray('e_z').GetTuple(q)
            Q = np.array((e_r,e_t,e_z))

            sigma_gnr_g = []
            stiffness_g = []
            varWallProps_g = []
            p_est_g = []

            inv_curr = self.vesselReference.GetCellData().GetArray('inv_curr').GetTuple1(q)

            for p in range(self.nG):
                sigma_gnr_mem = self.vesselReference.GetCellData().GetArray('sigma_mem').GetTuple(q)[p*9:(p+1)*9]
                sigma_gnr = rotateStress(sigma_gnr_mem,Q)
                sigma_gnr_g = sigma_gnr_g + sigma_gnr.tolist()
                p_est = (10.0*inv_curr-(sigma_gnr[0]+sigma_gnr[1]+sigma_gnr[2]))/3.0
                p_est_g = p_est_g + [p_est]

                J_di = self.vesselReference.GetCellData().GetArray('J_di').GetTuple(q)[p]
                stiffness_mem = self.vesselReference.GetCellData().GetArray('stiffness_mem').GetTuple(q)[p*36:(p+1)*36]
                stiffness = rotateStiffness(stiffness_mem,Q)
                stiffness_g = stiffness_g + stiffness.tolist()

                varWallProps_g = np.hstack((varWallProps_g,np.hstack((stiffness,p_est,J_di,sigma_gnr,p))))

            self.vesselReference.GetCellData().GetArray('varWallProps').SetTuple(q,varWallProps_g)
            self.vesselReference.GetCellData().GetArray('p_est').SetTuple(q,p_est_g)
            self.vesselReference.GetCellData().GetArray('sigma_gnr').SetTuple(q,sigma_gnr_g)

        return



    def initializeCylinder(self):
        print("Initializing vessel...")

        with open('FolderVesselConfigurationFiles/Native_in_handshake_') as f:
            lines = f.readlines()[1:]
            nativeIn = []
            for line in lines:
                nativeIn.append([float(x) for x in line.split()])

        #Build material array (constant for now)
        materialArray = [ nativeIn[7][0],nativeIn[4][0],nativeIn[4][1],nativeIn[4][2],nativeIn[2][0]*10.0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[5][2],nativeIn[2][4]*10.0,nativeIn[2][5],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][3],nativeIn[5][3],nativeIn[2][6]*10.0,nativeIn[2][7],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][4],nativeIn[5][4],nativeIn[2][8]*10.0,nativeIn[2][9],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][5],nativeIn[5][5],nativeIn[2][10]*10.0,nativeIn[2][11],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][1],nativeIn[5][1],nativeIn[2][2]*10.0,nativeIn[2][3],0,0,0,0,0,0,0,0,0,0,\
                          nativeIn[7][2],nativeIn[16][0],nativeIn[14][2],nativeIn[14][1],0,0,0,0,0,0,0,0,0,0]

        #********************

        points = np.empty([(self.numCirc+1)*(self.numLen+1)*(self.numRad+1),3])

        structureIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        innerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        outerIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        proximalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))
        distalIds = np.zeros((self.numCirc+1)*(self.numLen+1)*(self.numRad+1))

        proximalIds[0:(self.numCirc+1)*(self.numRad+1)]=1
        distalIds[(self.numCirc+1)*(self.numLen)*(self.numRad+1):(self.numCirc+1)*(self.numLen+1)*(self.numRad+1)]=1
        innerIds[0:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1
        outerIds[self.numRad:(self.numCirc+1)*(self.numLen+1)*(self.numRad+1):(self.numRad+1)]=1

        fluidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)
        solidStressQueryIds = np.zeros(self.numCirc*self.numLen*self.numRad)

        tevgValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        aneurysmValue = np.zeros(self.numCirc*self.numLen*self.numRad)
        simulateIds = np.ones(self.numCirc*self.numLen*self.numRad)

        e_r = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_t = np.zeros((self.numCirc*self.numLen*self.numRad,3))
        e_z = np.zeros((self.numCirc*self.numLen*self.numRad,3))

        e_ma = np.zeros((self.numCirc*self.numLen*self.numRad,98))

        for i in range(self.numLen+1):
            for j in range(self.numCirc+1):
                for k in range(self.numRad+1):
                    xPt = (self.radius + self.thickness*k/self.numRad)*np.cos(2.0*j*np.pi/self.numCirc)
                    yPt = (self.radius + self.thickness*k/self.numRad)*np.sin(2.0*j*np.pi/self.numCirc)
                    #zPt = -(1 - 2.0*i/(self.numLen))/abs(1 - 2.0*i/(self.numLen))*length*0.5*(np.exp(2.0*abs(1 - 2.0*i/(self.numLen)))-1.0)/(np.exp(2.0)-1.0)
                    zPt = self.length*i/self.numLen - self.length/2.0
                    structureIds[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k] = i*(self.numCirc)*(self.numRad+1) + (j%self.numCirc)*(self.numRad+1) + k + 1
                    points[i*(self.numCirc+1)*(self.numRad+1) + j*(self.numRad+1) + k,:] = [xPt, yPt, zPt]


        for i in range(self.numLen):
            for j in range(self.numCirc):
                for k in range(self.numRad):

                    xPt = np.cos(2.0*(j+0.5)*np.pi/self.numCirc)
                    yPt = np.sin(2.0*(j+0.5)*np.pi/self.numCirc)
                    zPt = self.length*(i+0.5)/self.numLen - self.length/2.0

                    v_r = np.array([xPt,yPt,0])
                    v_t = np.array([-yPt,xPt,0])
                    v_z = np.array([0,0,1])

                    e_r[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_r
                    e_t[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_t
                    e_z[i*self.numCirc*self.numRad + j*self.numRad + k,:] = v_z

                    A = np.array([v_r,v_t,v_z])
                    # Is this a rotation matrix?
                    if np.sometrue(np.abs(np.dot(np.array(A), np.transpose(np.array(A))) - 
                                          np.eye(3, dtype=float)) > 1e-6):
                        print(A)
                        raise RuntimeError('Matrix *A* does not describe a rotation.')

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,:] = materialArray

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,5:8]   = v_r
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,8:11]  = v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,11:14] = v_z

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,18:21] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,32:35] =  np.cos(0.0*np.pi/180.0)*v_z + np.sin(0.0*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,46:49] =  np.cos(41.94*np.pi/180.0)*v_z + np.sin(41.94*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,60:63] =  np.cos(318.06*np.pi/180.0)*v_z + np.sin(318.06*np.pi/180.0)*v_t

                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,74:77] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t
                    e_ma[k + j*self.numRad + i*self.numRad*self.numCirc,88:91] =  np.cos(90.0*np.pi/180.0)*v_z + np.sin(90.0*np.pi/180.0)*v_t

                    fluidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc
                    solidStressQueryIds[k + j*self.numRad + i*self.numRad*self.numCirc] = j*self.numRad + i*self.numRad*self.numCirc + self.numRad

                    if self.tevg:
                        tevgValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getTEVGValue([xPt,yPt,zPt], self.radius)
                    if self.aneurysm:
                        aneurysmValue[k + j*self.numRad + i*self.numRad*self.numCirc] = getAneurysmValue([xPt,yPt,zPt], self.radius)

        vol = pv.StructuredGrid()
        vol.points = points
        vol.dimensions = [self.numRad+1, self.numCirc+1, self.numLen+1]

        numCells = vol.GetNumberOfCells()
        vol.GetPointData().AddArray(pv.convert_array((proximalIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((distalIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((innerIds).astype(int),name="InnerRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outerIds).astype(int),name="OuterRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((structureIds).astype(int),name="StructureID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="GlobalElementID"))
        vol.GetCellData().AddArray(pv.convert_array(np.linspace(1,numCells,numCells).astype(int),name="CellStructureID"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="varWallProps"))

        #Exchange these to gauss point quantities
        #Matrices/arrays
        vol.GetCellData().AddArray(pv.convert_array(np.zeros((self.numCirc*self.numLen*self.numRad,45*8)).astype(float),name="temp_array"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.array([1,0,0,0,1,0,0,0,1]),8),(numCells,1)).astype(float),name="defGrad_mem"))

        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(36),8),(numCells,1)).astype(float),name="stiffness_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(9),8),(numCells,1)).astype(float),name="sigma_mem"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile(np.zeros(6),8),(numCells,1)).astype(float),name="sigma_gnr"))
        #Scalar values
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_curr"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([1],8),(numCells,1)).astype(float),name="J_di"))
        vol.GetCellData().AddArray(pv.convert_array(np.tile(np.tile([0],8),(numCells,1)).astype(float),name="p_est"))

        vol.GetCellData().AddArray(pv.convert_array(np.ones(numCells).astype(float),name="J_e"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][1]*np.ones(numCells).astype(float),name="wss_prev"))

        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_curr"))
        vol.GetCellData().AddArray(pv.convert_array(nativeIn[13][0]*np.ones(numCells).astype(float),name="inv_prev"))


        vol.GetCellData().AddArray(pv.convert_array(e_r.astype(float),name="e_r"))
        vol.GetCellData().AddArray(pv.convert_array(e_t.astype(float),name="e_t"))
        vol.GetCellData().AddArray(pv.convert_array(e_z.astype(float),name="e_z"))
        vol.GetCellData().AddArray(pv.convert_array(e_ma.astype(float),name="material_global"))
        vol.GetCellData().AddArray(pv.convert_array(fluidStressQueryIds.astype(int),name="fluidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(solidStressQueryIds.astype(int),name="solidStressQueryID"))
        vol.GetCellData().AddArray(pv.convert_array(simulateIds.astype(int),name="Simulate"))
        vol.GetCellData().AddArray(pv.convert_array(aneurysmValue.astype(float),name="aneurysmValue"))
        vol.GetCellData().AddArray(pv.convert_array(tevgValue.astype(float),name="tevgValue"))
        vol = vol.cast_to_unstructured_grid()
        vol = clean(vol, tolerance=1e-10)
        numPts = vol.GetNumberOfPoints()
        vol.GetPointData().AddArray(pv.convert_array(np.linspace(1,numPts,numPts).astype(int),name="GlobalNodeID"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="displacements"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_curr"))
        vol.GetPointData().AddArray(pv.convert_array(np.tile(np.zeros(3),(numPts,1)).astype(float),name="residual_prev"))
        vol.GetPointData().AddArray(pv.convert_array(np.zeros(numPts).astype(float),name="Pressure"))

        return vol



    def generateFluidMesh(self, surf):
        print("Making inner volume...")
        pid = 0
        numQuad = self.numCirc // 4 + 1
        numPerim = self.numCirc // 4
        quadFrac = 0.1
        numTrans = 40
        points = []

        inIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))
        outIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))
        wallIds = np.zeros((numQuad**2 + (numTrans-1)*self.numCirc)*(self.numLen+1))

        outIds[0:numQuad**2] = 1
        outIds[((self.numLen+1)*numQuad**2):((self.numLen+1)*numQuad**2)+ (numTrans-1)*self.numCirc] = 1

        inIds[((self.numLen+1)-1)*numQuad**2:(self.numLen+1)*numQuad**2] = 1
        inIds[-(numTrans-1)*self.numCirc:] = 1

        for ia in range(self.numLen+1):
            wallIds[((self.numLen+1)*numQuad**2) + (numTrans-2)*self.numCirc + ia*(numTrans-1)*(self.numCirc):((self.numLen+1)*numQuad**2) + (numTrans-1)*self.numCirc + ia*(numTrans-1)*(self.numCirc)] = 1


        numPoints = surf.GetNumberOfPoints()
        origIdsDict = dict(zip((surf.get_array('StructureID')-1).tolist(),range(numPoints)))

        # generate quadratic mesh
        for ia in range((self.numLen+1)):
            c1 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*0*(self.numCirc // 4)]))
            c2 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*1*(self.numCirc // 4)]))
            c3 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*2*(self.numCirc // 4)]))
            c4 = np.array(surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + (self.numRad+1)*3*(self.numCirc // 4)]))
            center = (c1 + c2 + c3 + c4)/4.0

            M = np.linalg.lstsq(np.array([[-0.5,-0.5,0,1],[0.5,-0.5,0,1],[0.5,0.5,0,1],[-0.5,0.5,0,1]]),np.vstack((c1,c2,c3,c4)))[0]
            for iy in range(numQuad):
                for ix in range(numQuad):

                    rady = quadFrac * ( iy / (numQuad - 1) - 0.5)
                    radx = quadFrac * ( ix / (numQuad - 1) - 0.5)

                    points.append(np.matmul([radx, rady, 0,1],M))
                    pid += 1

        # generate transition mesh
        for ia in range(self.numLen+1):
            for ir in range(1, numTrans):
                for ic in range(self.numCirc):
                    if ic <= numPerim:
                        quadPoint = points[ic + ia*(numQuad**2)]
                    elif numPerim < ic <= numPerim*2:
                        quadPoint = points[(ic-numPerim + 1)*numQuad - 1 + ia*(numQuad**2)]
                    elif numPerim*2 < ic <= numPerim*3:
                        quadPoint = points[numQuad**2 - 1 - (ic-numPerim*2) + ia*(numQuad**2)]
                    elif numPerim*3 < ic <= numPerim*4:
                        quadPoint = points[(numPerim-(ic-numPerim*3))*numQuad  + ia*(numQuad**2)]

                    outerPoint = surf.GetPoint(origIdsDict[(ia*(self.numRad+1)*self.numCirc) + ic*(self.numRad+1)])

                    points.append(((outerPoint - quadPoint) * ir/(numTrans - 1)) + quadPoint)
                    pid += 1


        coords=[[0, 1, 0],
                [0, 1, 1],
                [1, 1, 1],
                [1, 1, 0],
                [0, 0, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 0, 0]]

        cid = 0
        cells = []

        # generate quadratic mesh
        for ia in range(self.numLen):
            for iy in range(numQuad - 1):
                for ix in range(numQuad - 1):
                    ids = []
                    for c in coords:
                        ids += [(iy + c[0]) * numQuad + ix + c[1] + (ia + c[2]) * numQuad ** 2]
                    cells.append(ids)
                    cid += 1
     
        # generate transition mesh
        for ia in range(self.numLen):
            for ic in range(self.numCirc):
                ids = []
                for c in coords:
                        if c[1] == 1:
                            # circular side
                            ids += [(ic + c[0])%self.numCirc + (self.numLen+1) * numQuad ** 2 + (ia + c[2]) * (numTrans - 1) * self.numCirc]
                        else:
                            if (ic + c[0]) <= numPerim:
                                ids += [(ic + c[0]) + (ia + c[2])*(numQuad**2)]
                            elif numPerim < (ic + c[0]) <= numPerim*2:
                                ids += [((ic + c[0])-numPerim + 1)*numQuad - 1 + (ia + c[2])*(numQuad**2)]
                            elif numPerim*2 < (ic + c[0]) <= numPerim*3:
                                ids += [numQuad**2 - 1 - ((ic + c[0])-numPerim*2) + (ia + c[2])*(numQuad**2)]
                            elif numPerim*3 < (ic + c[0]) <= numPerim*4:
                                ids += [(numPerim-((ic + c[0])-numPerim*3))*numQuad  + (ia + c[2])*(numQuad**2)]
                cells.append(ids)
                cid += 1

        # generate circular mesh
        for ia in range(self.numLen):
            for ir in range(numTrans-2):
                for ic in range(self.numCirc):
                    ids = []
                    for c in coords:
                            ids += [(ic + c[0])%self.numCirc + (self.numLen+1) * numQuad ** 2 + (ia + c[2]) * (numTrans - 1) * self.numCirc + (ir + c[1]) * self.numCirc]
                    cells.append(ids)
                    cid += 1


        # each cell is a VTK_HEXAHEDRON
        cellArray = [ [8] + cell for cell in cells]
        cellTypes = np.empty(np.shape(cells)[0], dtype=np.uint8)
        cellTypes[:] = vtk.VTK_HEXAHEDRON
        cellOffset = range(0,np.shape(cells)[0]*9,9)

        vol = pv.UnstructuredGrid(np.array(cellOffset), np.array(cellArray), np.array(cellTypes), np.array(points))

        vol.GetPointData().AddArray(pv.convert_array((inIds).astype(int),name="ProximalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((outIds).astype(int),name="DistalRegionID"))
        vol.GetPointData().AddArray(pv.convert_array((wallIds).astype(int),name="OuterRegionID"))

        return vol



    """

    def setFlow()

    def setPressure()

    def setTimestep()

    def setDampening()

    def setPenalty()

    """