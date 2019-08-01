import os
import vtk, qt, ctk, slicer
import logging
from SegmentEditorEffects import *
import numpy as np
import math
import vtkSegmentationCorePython

class SegmentEditorEffect(AbstractScriptedSegmentEditorEffect):
  """This effect uses Watershed algorithm to partition the input volume"""

  def __init__(self, scriptedEffect):
    scriptedEffect.name = 'SRS-Filter'
    scriptedEffect.perSegment = True # this effect operates on all segments at once (not on a single selected segment)
    AbstractScriptedSegmentEditorEffect.__init__(self, scriptedEffect)

    self.logic = SRSFilterLogic(scriptedEffect)
    self.logic.logCallback = self.addLog

    # parameters
    self.parameters = []
    self.parameters.append({'max':50.0, 'min':0.0,'default':15.0, 'step':0.1, 'suffix':'mm', 'name':'Offset First Shrinkwrap:', 'id':ARG_OFFSETFIRSTSHRINKWRAP,'tooltip':''})
    self.parameters.append({'max':50.0, 'min':0.1,'default':10.0, 'step':0.1, 'suffix':'mm^3', 'name':'Spacing First Remesh:', 'id':ARG_SPACINGFIRSTREMESH,'tooltip':''})
    self.parameters.append({'max':10, 'min':0.1,'default':1.0, 'step':0.1, 'suffix':'mm^3', 'name':'Spacing Second Shrinkwrap:', 'id':ARG_SPACINGSECONDREMESH,'tooltip':''})
    self.parameters.append({'max':10, 'min':1,'default':3, 'step':1, 'suffix':'', 'name':'Iterations First Shrinkwrap:', 'id':ARG_ITERATIONSFIRSTSHRINKWRAP,'tooltip':''})
    self.parameters.append({'max':10, 'min':0,'default':5, 'step':1, 'suffix':'', 'name':'Iterations Second Shrinkwrap:', 'id':ARG_ITERATIONSSECONDSHRINKWRAP,'tooltip':''})
    self.parameters.append({'max':100.0, 'min':0.1,'default':20.0, 'step':0.1, 'suffix':'mm', 'name':'Raycast Search Edge Length:', 'id':ARG_RAYCASTSEARCHEDGELENGTH,'tooltip':''})
    self.parameters.append({'max':100.0, 'min':0.1,'default':2.0, 'step':0.1, 'suffix':'mm', 'name':'Raycast Output Edge Length:', 'id':ARG_RAYCASTOUTPUTEDGELENGTH,'tooltip':''})
    self.parameters.append({'max':50.0, 'min':0.1,'default':2.0, 'step':0.01, 'suffix':'mm', 'name':'Raycast Max. Hit Distance:', 'id':ARG_RAYCASTMAXHITDISTANCE,'tooltip':''})
    self.parameters.append({'max':1000.0, 'min':0.1,'default':100.0, 'step':0.1, 'suffix':'mm', 'name':'Raycast Max. Length:', 'id':ARG_RAYCASTMAXLENGTH,'tooltip':''})
    self.parameters.append({'max':1000.0, 'min':0.0,'default':0.0, 'step':0.1, 'suffix':'mm', 'name':'Raycast Min. Length:', 'id':ARG_RAYCASTMINLENGTH,'tooltip':''})
    self.parameters.append({'max':10, 'min':0.01,'default':0.7, 'step':0.01, 'suffix':'mm', 'name':'Max. Models Distance:', 'id':ARG_MAXMODELSDISTANCE,'tooltip':''})
    self.parameters.append({'max':20.0, 'min':0.1,'default':1.5, 'suffix':'mm', 'step':0.1, 'name':'Solidification Thickness:', 'id':ARG_SOLIDIFICATIONTHICKNESS,'tooltip':''})
    self.parameters.append({'max':1.0, 'min':0.0,'default':0.2, 'suffix':'', 'step':0.1, 'name':'Smoothing Factor:', 'id':ARG_SMOOTHINGFACTOR,'tooltip':''})
    
    # filter modes
    self.filterModes = []
    self.filterModes.append({'name':'Convex Hull', 'id':MODE_CONVEXHULL, 'default':False})
    self.filterModes.append({'name':'Raycasts', 'id':MODE_RAYCASTS,'default':False})
    self.filterModes.append({'name':'Deep Hull', 'id':MODE_DEEPHULL,'default':False})
    self.filterModes.append({'name':'Nonmanifold', 'id':MODE_NONMANIFOLD,'default':False})
    self.filterModes.append({'name':'Solidified Surface', 'id':MODE_SOLIDIFIED,'default':True})

    # outputTypes
    self.outputTypes = []
    self.outputTypes.append({'name':'Model', 'id':OUTPUT_MODEL, 'default':True})
    self.outputTypes.append({'name':'Segmentation', 'id':OUTPUT_SEGMENTATION, 'default':False})


  def clone(self):
    # It should not be necessary to modify this method
    import qSlicerSegmentationsEditorEffectsPythonQt as effects
    clonedEffect = effects.qSlicerSegmentEditorScriptedEffect(None)
    clonedEffect.setPythonSource(__file__.replace('\\','/'))
    return clonedEffect

  def icon(self):
    # It should not be necessary to modify this method
    iconPath = os.path.join(os.path.dirname(__file__), 'SegmentEditorEffect.png')
    if os.path.exists(iconPath):
      return qt.QIcon(iconPath)
    return qt.QIcon()

  def helpText(self):
    return """<html>This Filter results in an even surface of the current segment<br>.
    It is using a combination of shrinkwrapping, remeshing, raycasting and solidification algorithms. Because of the process, the following intermediate outputs are possible (with inheritting parameter dependencies):<br>

    <ol style="margin: 0">
    <li>Convex Hull: Offset First Shrinkwrap, Spacing First Remesh, Iterations First Shrinkwrap.</li>
    <li>Raycast Result (especially for parameter testing reasons): Raycast Search/Output Edge Lenght, Raycast Max. Hit Distance, Raycast Max./Min. Length.</li>
    <li>Deep Hull: Spacing Second Remesh, Iterations Second Shrinkwrap.</li>
    <li>Nonsolid Model (especially for further processing in CAD software): Max. Model Distance, Smoothing Factor. WARNING: It is not possible to create a segment out of this nonsolid process result.</li>
    <li>Solidified Surface: Solidification Thickness.</li>
    </ol><br>
    
    Current parameters are especially fitted for working on fractured hemipelvic bone segmentation. For further information, license, disclaimers and possible research partnerships visit <a href="https://github.com/sebastianandress/Slicer-SegmentEditorSRSFilter">this</a> github repository.
    </html>"""

  def activate(self):
    pass

  def deactivate(self):
    self.cleanup()
  
  def cleanup(self):
    pass


  def setupOptionsFrame(self):
    
    # setup filter modes

    self.filterModeLayout = qt.QVBoxLayout()
    self.scriptedEffect.addLabeledOptionsWidget("Filter Mode:", self.filterModeLayout)
    self.filterModeGroup = qt.QButtonGroup()

    for i, mode in enumerate(self.filterModes):
      mode['element'] = qt.QRadioButton(mode['name'])
      self.filterModeGroup.addButton(mode['element'])
      self.filterModeLayout.addWidget(mode['element'])
      if mode['default']: mode['element'].setChecked(True)
    
    self.filterModeGroup.connect('buttonClicked(int)', self.updateMRMLFromGUI)

    # setup output type
    self.outputTypeLayout = qt.QVBoxLayout()
    self.scriptedEffect.addLabeledOptionsWidget("Output Type:", self.outputTypeLayout)
    self.outputTypeGroup = qt.QButtonGroup()

    for i, out in enumerate(self.outputTypes):
      out['element'] = qt.QRadioButton(out['name'])
      self.outputTypeGroup.addButton(out['element'])
      self.outputTypeLayout.addWidget(out['element'])
      if out['default']: out['element'].setChecked(True)

    self.outputTypeGroup.connect('buttonClicked(int)', self.updateMRMLFromGUI)


    # setup parameters
    advancedSettingsFrame = ctk.ctkCollapsibleGroupBox()
    advancedSettingsFrame.objectName = 'advancedSettingsFrame'
    advancedSettingsFrame.title = "Advanced Settings"
    advancedSettingsFrame.collapsed = True
    advancedSettingsFrame.setLayout(qt.QFormLayout())

    for param in self.parameters:
      param['element'] = slicer.qMRMLSliderWidget()
      param['element'].setMRMLScene(slicer.mrmlScene)
      #param['element'].quantity = "length" # get unit, precision, etc. from MRML unit node
      param['element'].minimum = param['min']
      param['element'].maximum = param['max']
      param['element'].singleStep = param['step']
      param['element'].value = param['default']
      param['element'].suffix = param['suffix']
      param['element'].setToolTip(param['tooltip'])
      #self.scriptedEffect.addLabeledOptionsWidget(param['name'], param['element'])
      advancedSettingsFrame.layout().addRow(param['name'], param['element'])
      param['element'].connect('valueChanged(double)', self.updateMRMLFromGUI)
    
    self.scriptedEffect.addOptionsWidget(advancedSettingsFrame)

    # Apply Button

    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.objectName = self.__class__.__name__ + 'Apply'
    self.applyButton.setToolTip("")
    self.scriptedEffect.addOptionsWidget(self.applyButton)

    self.applyButton.connect('clicked()', self.onApply)


  def createCursor(self, widget):
    return slicer.util.mainWindow().cursor

  def layoutChanged(self):
    pass

  def processInteractionEvents(self, callerInteractor, eventId, viewWidget):
    return False # For the sake of example

  def processViewNodeEvents(self, callerViewNode, eventId, viewWidget):
    pass # For the sake of example

  def setMRMLDefaults(self):
    for param in self.parameters:
      self.scriptedEffect.setParameterDefault(param['id'], param['default'])

    self.scriptedEffect.setParameterDefault(ARG_FILTERMODE, next(item for item in self.filterModes if item['default'] == True)['id'])
    self.scriptedEffect.setParameterDefault(ARG_OUTPUTTYPE, next(item for item in self.outputTypes if item['default'] == True)['id'])

  def updateGUIFromMRML(self):
    for param in self.parameters:
      value = self.scriptedEffect.doubleParameter(param['id'])
      wasBlocked = param['element'].blockSignals(True)
      param['element'].value = abs(value)
      param['element'].blockSignals(wasBlocked)
    
    filterModeID = self.scriptedEffect.parameter(ARG_FILTERMODE)
    for l in self.filterModes:
      if l['id'] == filterModeID:
        l['element'].setChecked(True)
      else:
        l['element'].setChecked(False)
    
    outputTypeID = self.scriptedEffect.parameter(ARG_OUTPUTTYPE)
    
    for t in self.outputTypes:
      if t['id'] == outputTypeID:
        t['element'].setChecked(True)
      else:
        t['element'].setChecked(False)

    self.disableOptionsLayout()
    self.cleanup()
    

  def updateMRMLFromGUI(self):
    for param in self.parameters:
      self.scriptedEffect.setParameter(param['id'], param['element'].value)

    for mode in self.filterModes:
      if mode['element'].isChecked():
        self.scriptedEffect.setParameter(ARG_FILTERMODE, mode['id'])
    
    for typ in self.outputTypes:
      if typ['element'].isChecked():
        self.scriptedEffect.setParameter(ARG_OUTPUTTYPE, typ['id'])

    self.disableOptionsLayout()
    self.cleanup()
  

  #
  # Effect specific methods (the above ones are the API methods to override)
  #

  def disableOptionsLayout(self):
    [f['element'].setEnabled(True) for f in self.filterModes if f.has_key('element')]
    [f['element'].setEnabled(True) for f in self.outputTypes if f.has_key('element')]
    [f['element'].setEnabled(True) for f in self.parameters if f.has_key('element')]
    if self.scriptedEffect.parameter(ARG_FILTERMODE) == MODE_SOLIDIFIED: return
    
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_SOLIDIFICATIONTHICKNESS]
    if self.scriptedEffect.parameter(ARG_FILTERMODE) == MODE_NONMANIFOLD:
      [f['element'].setEnabled(False) for f in self.outputTypes if f['id'] == OUTPUT_SEGMENTATION]
      [f['element'].setChecked(True) for f in self.outputTypes if f['id'] == OUTPUT_MODEL]
      return
    
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_MAXMODELSDISTANCE]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_SMOOTHINGFACTOR]
    if self.scriptedEffect.parameter(ARG_FILTERMODE) == MODE_DEEPHULL: return
    
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_ITERATIONSSECONDSHRINKWRAP]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_SPACINGSECONDREMESH]    
    if self.scriptedEffect.parameter(ARG_FILTERMODE) == MODE_RAYCASTS: return

    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_RAYCASTMAXHITDISTANCE]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_RAYCASTMAXLENGTH]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_RAYCASTMINLENGTH]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_RAYCASTOUTPUTEDGELENGTH]
    [f['element'].setEnabled(False) for f in self.parameters if f['id'] == ARG_RAYCASTSEARCHEDGELENGTH]
    return


  def onApply(self):

    self.scriptedEffect.saveStateForUndo()
    
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

    #TODO: check if segment is not empty

    seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    segID = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()

    kwargs = {}
    for arg in self.parameters:
      kwargs.update({arg['id']:self.scriptedEffect.doubleParameter(arg['id'])})
    kwargs.update({ARG_FILTERMODE:self.scriptedEffect.parameter(ARG_FILTERMODE)})
    kwargs.update({ARG_OUTPUTTYPE:self.scriptedEffect.parameter(ARG_OUTPUTTYPE)})
    self.logic.ApplySRSFilter(seg, segID, **kwargs)

    qt.QApplication.restoreOverrideCursor()


  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update



class SRSFilterLogic(object):

  def __init__(self, scriptedEffect):
    self.scriptedEffect = scriptedEffect
    self.logCallback = None

  def ApplySRSFilter(self, segmentationNode, segmentID, **kwargs):
    """[summary]
    
    Arguments:
        segmentationNode {[vtkMRMLSegmentationNode]} -- [description]
        segmentID {[string]} -- [description]
        **outputType {[string]} -- [description]
        **filterMode {[string]} -- [description]
        **offsetFirstShrinkwrap {[double]} -- []
        **spacingFirstRemesh {[double]} -- []
        **spacingSecondRemesh {[double]} -- []
        **iterationsFirstShrinkwrap {[int]} -- []
        **iterationsSecondShrinkwrap {[int]} -- []
        **raycastSearchEdgeLength {[double]} -- []
        **raycastOutputEdgeLength {[double]} -- []
        **raycastMaxHitDistance {[double]} -- []
        **raycastMaxLength {[double]} -- []
        **raycastMinLength {[double]} -- []
        **maxModelDistance {[double]} -- []
        **solidificationThickness {[double]} -- []
        **smoothingFactor {[double]} -- []

    
    Returns:
        [bool] -- [description]
    """

    
    options = {
      ARG_OUTPUTTYPE : OUTPUT_MODEL,
      ARG_FILTERMODE : MODE_SOLIDIFIED,
      ARG_OFFSETFIRSTSHRINKWRAP : 15,
      ARG_SPACINGFIRSTREMESH : 10,
      ARG_SPACINGSECONDREMESH : 1,
      ARG_ITERATIONSFIRSTSHRINKWRAP : 3,
      ARG_ITERATIONSSECONDSHRINKWRAP : 5,
      ARG_RAYCASTSEARCHEDGELENGTH : 20,
      ARG_RAYCASTOUTPUTEDGELENGTH : 2,
      ARG_RAYCASTMAXHITDISTANCE : 2,
      ARG_RAYCASTMAXLENGTH : 100,
      ARG_RAYCASTMINLENGTH : 0,
      ARG_MAXMODELSDISTANCE : 0.7,
      ARG_SOLIDIFICATIONTHICKNESS : 1.5,
      ARG_SMOOTHINGFACTOR : 0.2
      }

    options.update(kwargs)

    self.segLogic = slicer.vtkSlicerSegmentationsModuleLogic
    self.modelsLogic = slicer.modules.models.logic()

    tempNodes = []

    def polydataToModel(polydata, smooth=True):
      if self.logCallback: self.logCallback('Creating Model...')
      
      if smooth:
        polydata = smoothPolydata(polydata)

      modelNode = self.modelsLogic.AddModel(polydata)
        
      seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode().GetSegmentation().GetSegment(self.scriptedEffect.parameterSetNode().GetSelectedSegmentID())
      modelNode.GetDisplayNode().SetColor(seg.GetColor())
      modelNode.SetName(seg.GetName())
      modelNode.GetDisplayNode().SliceIntersectionVisibilityOn()

      return True

    def polydataToSegment(polydata, smooth=True):
      if self.logCallback: self.logCallback('Updating Segmentation...')
      
      if smooth:
        polydata = smoothPolydata(polydata)

      tempSegment = vtkSegmentationCorePython.vtkSegment()
      tempSegment.SetName(segmentationNode.GetSegmentation().GetSegment(segmentID).GetName())
      tempSegment.SetColor(segmentationNode.GetSegmentation().GetSegment(segmentID).GetColor())
      tempSegment.AddRepresentation(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(), polydata)
      segmentationNode.GetSegmentation().GetSegment(segmentID).DeepCopy(tempSegment)
      
      # Update views
      segmentationNode.Modified()
      if segmentationNode.GetDisplayNode():
        segmentationNode.GetDisplayNode().SetPreferredDisplayRepresentationName3D(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName())
      segmentationNode.CreateClosedSurfaceRepresentation()

      return True

    def cleanup():
      for node in tempNodes:
        if node:
          slicer.mrmlScene.RemoveNode(node)
      
      if self.logCallback: self.logCallback('')

    def remeshPolydata(polydata, spacing):
      whiteImage = vtk.vtkImageData()
      bounds = [0]*6
      polydata.GetBounds(bounds)

      #spacing = [spacing]*3
      whiteImage.SetSpacing(spacing)

      dim = [0]*3
      for i in range(3):
        dim[i] = int(math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) + 1
        if (dim[i] < 1):
          dim[i] = 1
      whiteImage.SetDimensions(dim)
      #whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)
      whiteImage.SetExtent(0, dim[0], 0, dim[1], 0, dim[2])
      origin = [0]*3

      origin[0] = bounds[0]# + spacing[0] / 2
      origin[1] = bounds[2]# + spacing[1] / 2
      origin[2] = bounds[4]# + spacing[2] / 2
      whiteImage.SetOrigin(origin)

      whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)

      pol2stenc = vtk.vtkPolyDataToImageStencil()
      pol2stenc.SetInputData(polydata)

      pol2stenc.SetOutputOrigin(origin)
      pol2stenc.SetOutputSpacing(spacing)
      pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
      pol2stenc.Update()

      imgstenc = vtk.vtkImageStencil()
      imgstenc.SetInputData(whiteImage)
      imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
      imgstenc.ReverseStencilOn()
      imgstenc.SetBackgroundValue(1)
      imgstenc.Update()

      revimgstenc = vtk.vtkImageStencil()
      revimgstenc.SetInputData(imgstenc.GetOutput())
      revimgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
      revimgstenc.ReverseStencilOff()
      revimgstenc.SetBackgroundValue(0)
      revimgstenc.Update()

      discreteCubes = vtk.vtkDiscreteMarchingCubes()
      discreteCubes.SetInputConnection(revimgstenc.GetOutputPort())
      discreteCubes.GenerateValues(1,0,0)
      discreteCubes.Update()

      reverse = vtk.vtkReverseSense()
      reverse.SetInputConnection(discreteCubes.GetOutputPort())
      reverse.ReverseCellsOn()
      reverse.ReverseNormalsOn()
      reverse.Update()

      return reverse.GetOutput()

    def smoothPolydata(polydata):
      passBand = pow(10.0, -4.0*options[ARG_SMOOTHINGFACTOR])
      smootherSinc = vtk.vtkWindowedSincPolyDataFilter()
      # smootherSinc.SetPassBand(passBand)
      # smootherSinc.SetInputConnection(decimator.GetOutputPort())
      smootherSinc.SetInputData(polydata)
      smootherSinc.SetNumberOfIterations(20)
      smootherSinc.FeatureEdgeSmoothingOff()
      smootherSinc.BoundarySmoothingOff()
      # smootherSinc.ReleaseDataFlagOn()
      smootherSinc.NonManifoldSmoothingOn()
      smootherSinc.NormalizeCoordinatesOn()
      smootherSinc.Update()

      return smootherSinc.GetOutput()


    if self.logCallback: self.logCallback('Filtering process started...')

    segmentationNode.CreateClosedSurfaceRepresentation()
    segmentationNode.CreateBinaryLabelmapRepresentation()
    segment = segmentationNode.GetSegmentation().GetSegment(segmentID)
    inputPolyData = segment.GetRepresentation(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName())
    inputSegmentLabelmap = segment.GetRepresentation(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationBinaryLabelmapRepresentationName())


    #region create sphere
    bounds = np.array([0]*6)
    inputPolyData.GetBounds(bounds)
    dimensions = np.array([bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])
    center = np.array([bounds[0]+dimensions[0]/2, bounds[2]+dimensions[1]/2, bounds[4]+dimensions[2]/2])
    radius = max(dimensions)

    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetRadius(radius)
    sphereSource.SetPhiResolution(int(radius/10))
    sphereSource.SetThetaResolution(int(radius/10))
    sphereSource.SetCenter(center)
    sphereSource.Update()

    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputConnection(sphereSource.GetOutputPort())
    cleanPolyData.Update()

    shrinkModelPD = vtk.vtkPolyData()
    shrinkModelPD.DeepCopy(cleanPolyData.GetOutput())


    #endregion

    #region Shrinkwrap

    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(inputPolyData)
    cellLocator.BuildLocator()
  
    for x in range(int(options[ARG_ITERATIONSFIRSTSHRINKWRAP])):
      
      # shrinkwrap
      if self.logCallback: self.logCallback('Shrinkwrapping %s/%s...' %(x+1, int(options[ARG_ITERATIONSFIRSTSHRINKWRAP])))
      
      points = shrinkModelPD.GetPoints()

      for i in xrange(points.GetNumberOfPoints()):
        originPoint = np.array(points.GetPoint(i))
        closestPoint = np.array([0.0,0.0,0.0])
        cell = vtk.vtkGenericCell()
        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        closestPointDist2 = vtk.mutable(0)

        cellLocator.FindClosestPoint(originPoint, closestPoint, cell, cellId, subId, closestPointDist2)

        vector = closestPoint - originPoint
        vectorLength = np.linalg.norm(vector)

        if options[ARG_OFFSETFIRSTSHRINKWRAP] > 0 and vectorLength > 0.01:
          newLocation = closestPoint - ((vector/vectorLength) * options[ARG_OFFSETFIRSTSHRINKWRAP])
        else:
          newLocation = closestPoint
        
        points.SetPoint(i, newLocation)
      
      shrinkModelPD.SetPoints(points)

      if x == (int(options[ARG_ITERATIONSFIRSTSHRINKWRAP]) - 1):
        break

      # remesh
      if self.logCallback: self.logCallback('Remeshing %s/%s...' %(x+1, int(options[ARG_ITERATIONSSECONDSHRINKWRAP])))
      shrinkModelPD.DeepCopy(remeshPolydata(shrinkModelPD, [options[ARG_SPACINGFIRSTREMESH]]*3))
    
    if options[ARG_FILTERMODE] == MODE_CONVEXHULL:
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        polydataToSegment(shrinkModelPD, False)
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        polydataToModel(shrinkModelPD, False)
      else:
        logging.error('unknown outputType')
        return False
      
      cleanup()
      return True

    #endregion

    #region Raycast
    if self.logCallback: self.logCallback('Raycasting...')

    # Find Large Faces and remember IDs of connected points
    largeCellIds = vtk.vtkIdList() # IDs of cells
    for i in xrange(shrinkModelPD.GetNumberOfCells()):
      cell = shrinkModelPD.GetCell(i)

      # get Length longest edge of cell
      pointsArray = list()
      for p in xrange(cell.GetNumberOfPoints()):
        pointsArray.append(np.array(cell.GetPoints().GetPoint(p)))

      edgeLength = list()
      for pa in xrange(len(pointsArray) - 1):
        length = np.linalg.norm(pointsArray[pa] - pointsArray[pa + 1])
        edgeLength.append(length)

      if max(edgeLength) > options[ARG_RAYCASTSEARCHEDGELENGTH]:
        largeCellIds.InsertNextId(i)

    # extract large cells for cell point localization

    largeCellsPolyData = vtk.vtkPolyData()
    largeCellsPolyData.DeepCopy(shrinkModelPD)
    largeCellsPolyData.BuildLinks()

    for c in xrange(largeCellsPolyData.GetNumberOfCells()):
      if largeCellIds.IsId(c) == -1:
        largeCellsPolyData.DeleteCell(c)

    largeCellsPolyData.RemoveDeletedCells()

    # subdivide
    ids = vtk.vtkIdFilter()
    adapt = vtk.vtkAdaptiveSubdivisionFilter()
    adapt.SetInputData(shrinkModelPD)
    adapt.SetMaximumEdgeLength(options[ARG_RAYCASTOUTPUTEDGELENGTH])
    adapt.SetMaximumTriangleArea(vtk.VTK_INT_MAX)
    adapt.SetMaximumNumberOfPasses(vtk.VTK_INT_MAX)
    adapt.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(adapt.GetOutput())
    clean.Update()

    shrinkModelPD.DeepCopy(clean.GetOutput())

    if largeCellIds.GetNumberOfIds() > 0 and options[ARG_RAYCASTMAXLENGTH] > 0.0:

      # locate the points of previous large cells and write into largePointIds Set
      largeDistance = vtk.vtkImplicitPolyDataDistance()
      largeDistance.SetInput(largeCellsPolyData)

      largePointIds = set()
      for p in xrange(shrinkModelPD.GetNumberOfPoints()):
        point = [0.0]*3
        shrinkModelPD.GetPoints().GetPoint(p, point)
        distance = largeDistance.EvaluateFunction(point)
        if abs(distance) < 1:
          largePointIds.update([p])

      # generate normals

      normals = vtk.vtkPolyDataNormals()
      normals.ComputePointNormalsOn()
      normals.ComputeCellNormalsOff()
      normals.SplittingOff()
      normals.SetInputData(shrinkModelPD)
      normals.AutoOrientNormalsOff()
      normals.FlipNormalsOff()
      normals.Update()

      shrinkModelPD.DeepCopy(normals.GetOutput())
    
      # projection
      vert_location_dict = {} # dict to save all projection results

      for i in xrange(shrinkModelPD.GetNumberOfCells()):
        cell = shrinkModelPD.GetCell(i)
        pointIds = cell.GetPointIds()
        
        if cell.GetCellType() == 5: # cell with face
          for p in xrange(pointIds.GetNumberOfIds()):
            pointId = pointIds.GetId(p)

            # check if cell with point was large before subdividion, and if point got checked already
            if not pointId in largePointIds or (vert_location_dict.has_key(pointId) and vert_location_dict[pointId][0] == False):
              
              # random False value, won't be moved
              vert_location_dict.update({pointId:(False,np.array([0.0,0.0,0.0]),0.0)})

            else:
              cell = shrinkModelPD.GetCell(i)
              pointId = cell.GetPointIds().GetId(p)
              normal = np.array(shrinkModelPD.GetPointData().GetArray('Normals').GetTuple(pointId)) * (-1)
              vector = normal * options[ARG_RAYCASTMAXLENGTH] # max Length greater 0, checked above

              points = cell.GetPoints()

              # find intersection (point + vector) and label model
              a0 = points.GetPoint(p)
              a1 = a0 + vector
              tol = 1.0

              t = vtk.mutable(0)
              glo = np.array([0.0,0.0,0.0]) #global
              par = np.array([0.0,0.0,0.0]) #parametric
              cell = vtk.vtkGenericCell()
              cellId = vtk.mutable(0)
              subId = vtk.mutable(0)
              cellLocator.IntersectWithLine(a0, a1, tol, t, glo, par, subId, cellId, cell)

              loc_new = np.array(glo)# - (normal * options[ARG_OFFSETFIRSTSHRINKWRAP])
              length = np.linalg.norm(glo - a0)
              res = False
              if np.linalg.norm(glo) != 0:
                res = True
              vert_location_dict.update({pointId:(res,loc_new,length)})

      numberOfPoints = shrinkModelPD.GetNumberOfPoints()

      for i in xrange(numberOfPoints):
        # check result
        if vert_location_dict[i][0] == True:
          # check min length
          if vert_location_dict[i][2] > options[ARG_RAYCASTMINLENGTH]:
            
            # check distance between two new locations with positive result
            cellIds = vtk.vtkIdList()
            shrinkModelPD.GetPointCells(i, cellIds)
            pointChanged = False
            for c in xrange(cellIds.GetNumberOfIds()):
              if pointChanged == True:
                break
              cell = shrinkModelPD.GetCell(cellIds.GetId(c))
              pointIds = cell.GetPointIds()
              for p in xrange(pointIds.GetNumberOfIds()):
                if pointChanged == True:
                  break
                pointId = pointIds.GetId(p)
                if pointId != i and vert_location_dict[pointId][0] == True:
                  point = vert_location_dict[pointId][1]
                  distance = np.linalg.norm(-vert_location_dict[i][1] + point)
                  if distance < options[ARG_RAYCASTMAXHITDISTANCE]:
                    shrinkModelPD.GetPoints().SetPoint(i, vert_location_dict[i][1])
                    pointChanged = True

    if options[ARG_FILTERMODE] == MODE_RAYCASTS:
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        polydataToSegment(shrinkModelPD, False)
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        polydataToModel(shrinkModelPD, False)
      else:
        logging.error('unknown outputType')
        return False
      
      cleanup()
      return True
    
    #endregion

    #region Shrinkwrap

    for x in range(int(options[ARG_ITERATIONSSECONDSHRINKWRAP])+1):
      
      # remesh
      if self.logCallback: self.logCallback('Remeshing %s/%s...' %(x+1, int(options[ARG_ITERATIONSSECONDSHRINKWRAP])+1))
      shrinkModelPD.DeepCopy(remeshPolydata(shrinkModelPD, [options[ARG_SPACINGSECONDREMESH]]*3))

      if x == int(options[ARG_ITERATIONSSECONDSHRINKWRAP]):
        break

      # shrinkwrap
      if self.logCallback: self.logCallback('Shrinkwrapping %s/%s...' %(x+1, int(options[ARG_ITERATIONSSECONDSHRINKWRAP])))

      smoothFilter = vtk.vtkSmoothPolyDataFilter()
      smoothFilter.SetInputData(0, shrinkModelPD)
      smoothFilter.SetInputData(1, inputPolyData)
      smoothFilter.Update()
      shrinkModelPD.DeepCopy(smoothFilter.GetOutput())

    shrinkModelPD.DeepCopy(smoothPolydata(shrinkModelPD))
    
    if options[ARG_FILTERMODE] == MODE_DEEPHULL:
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        polydataToSegment(shrinkModelPD, False)
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        polydataToModel(shrinkModelPD, False)
      else:
        logging.error('unknown outputType')
        return False
      
      cleanup()
      return True

    #endregion

    #region Remove Caps
    if self.logCallback: self.logCallback('Removing Caps...')

    # implicit distance, add point ids with larger distance to ids
    implicitDistance = vtk.vtkImplicitPolyDataDistance()
    implicitDistance.SetInput(inputPolyData)
    
    # delete cells in great distance
    nonsolidPolyData = vtk.vtkPolyData()
    nonsolidPolyData.DeepCopy(shrinkModelPD)
    nonsolidPolyData.BuildLinks()

    for c in range(nonsolidPolyData.GetNumberOfCells()):
      cell = nonsolidPolyData.GetCell(c)
      points = cell.GetPoints()
      for p in xrange(points.GetNumberOfPoints()):
        point = points.GetPoint(p)
        distance = implicitDistance.EvaluateFunction(point)

        if abs(distance) > options[ARG_MAXMODELSDISTANCE]:
          nonsolidPolyData.DeleteCell(c)
          break

    nonsolidPolyData.RemoveDeletedCells()
    shrinkModelPD.DeepCopy(nonsolidPolyData)

    if options[ARG_FILTERMODE] == MODE_NONMANIFOLD:
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        logging.error('outputType "SEGMENTATION" not possible for filterMode "NONMANIFOLD".')
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        polydataToModel(shrinkModelPD, False)
      else:
        logging.error('unknown outputType')
        return False
      
      cleanup()
      return True

    #endregion

    #region Solidification
    if self.logCallback: self.logCallback('Solidifying...')

    # remove double vertices
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(shrinkModelPD)
    cleanPolyData.Update()

    # create normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetComputeCellNormals(1)
    normals.SetInputData(cleanPolyData.GetOutput())
    normals.SplittingOff()
    normals.Update()

    #shrinkModelPD = vtk.vtkPolyData()
    shrinkModelPD.DeepCopy(normals.GetOutput())
    numberOfPoints = shrinkModelPD.GetNumberOfPoints()

    # get boundary edges, used later
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.BoundaryEdgesOn()
    featureEdges.ColoringOff()
    featureEdges.FeatureEdgesOff()
    featureEdges.NonManifoldEdgesOff()
    featureEdges.ManifoldEdgesOff()
    featureEdges.SetInputData(normals.GetOutput())
    featureEdges.Update()

    addingPoints = []
    addingPolys = []


    for pointID in range(numberOfPoints):
      cellIDs = vtk.vtkIdList()
      shrinkModelPD.GetPointCells(pointID, cellIDs)
      normalsArray = []

      
      # ilterate through all cells/faces which contain point
      for id in range(cellIDs.GetNumberOfIds()):
        n = []
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3))
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3 + 1))
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3 + 2))

        normalsArray.append(np.array(n) * (-1))

      # calculate position of new vert
      dir_vec = np.zeros(3)
      
      for n in normalsArray:
        dir_vec = dir_vec + np.array(n)

      dir_vec_norm = dir_vec / np.linalg.norm(dir_vec)
      proj_length = np.dot(dir_vec_norm, np.array(normalsArray[0]))
      dir_vec_finallenght = dir_vec_norm * proj_length
      vertex_neu = np.array(shrinkModelPD.GetPoint(pointID)) + (dir_vec_finallenght * options[ARG_SOLIDIFICATIONTHICKNESS])
      
      # append point
      addingPoints.append(vertex_neu)

    for cellID in range(shrinkModelPD.GetNumberOfCells()):
      pointIDs = vtk.vtkIdList()
      shrinkModelPD.GetCellPoints(cellID, pointIDs)

      newPointIDs = vtk.vtkIdList()
      for i in reversed(range(pointIDs.GetNumberOfIds())):
        newPointIDs.InsertNextId(int(pointIDs.GetId(i) + numberOfPoints))

      addingPolys.append(newPointIDs)

    doubleSurfacePoints = vtk.vtkPoints()
    doubleSurfacePolys = vtk.vtkCellArray()

    doubleSurfacePoints.DeepCopy(shrinkModelPD.GetPoints())
    doubleSurfacePolys.DeepCopy(shrinkModelPD.GetPolys())

    for p in addingPoints:
      doubleSurfacePoints.InsertNextPoint(p)
    for p in addingPolys:
      doubleSurfacePolys.InsertNextCell(p)

    doubleSurfacePD = vtk.vtkPolyData()
    doubleSurfacePD.SetPoints(doubleSurfacePoints)
    doubleSurfacePD.SetPolys(doubleSurfacePolys)

    # add faces to boundary edges
    mergePoints = vtk.vtkMergePoints()
    mergePoints.InitPointInsertion(doubleSurfacePD.GetPoints(), doubleSurfacePD.GetBounds())
    mergePoints.SetDataSet(doubleSurfacePD)
    mergePoints.BuildLocator()

    manifoldPolys = vtk.vtkCellArray()
    manifoldPolys.DeepCopy(doubleSurfacePD.GetPolys())
    manifoldPoints = vtk.vtkPoints()
    manifoldPoints.DeepCopy(doubleSurfacePD.GetPoints())

    for e in range(featureEdges.GetOutput().GetNumberOfCells()):
      pointIDs = vtk.vtkIdList()
      featureEdges.GetOutput().GetCellPoints(e, pointIDs)
      if pointIDs.GetNumberOfIds() == 2: # -> Edge
        matchingPointIDs = []
        newPointIDs = vtk.vtkIdList()
        for p in range(2):
          matchingPointIDs.append(mergePoints.IsInsertedPoint(featureEdges.GetOutput().GetPoint(pointIDs.GetId(p))))
        if not (-1) in matchingPointIDs: # edge vertex not found in original pd, should not happen
          newPointIDs.InsertNextId(matchingPointIDs[1])
          newPointIDs.InsertNextId(matchingPointIDs[0])
          newPointIDs.InsertNextId(matchingPointIDs[0]+numberOfPoints)
          newPointIDs.InsertNextId(matchingPointIDs[1]+numberOfPoints)
          manifoldPolys.InsertNextCell(newPointIDs)

    manifoldPD = vtk.vtkPolyData()
    manifoldPD.SetPoints(manifoldPoints)
    manifoldPD.SetPolys(manifoldPolys)

    triangleFilter = vtk.vtkTriangleFilter()
    triangleFilter.SetInputData(manifoldPD)
    triangleFilter.Update()

    #endregion

    if options[ARG_FILTERMODE] == MODE_SOLIDIFIED:
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        polydataToSegment(triangleFilter.GetOutput(), False)
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        polydataToModel(triangleFilter.GetOutput(), False)
      else:
        logging.error('unknown outputType')
        return False
      
      cleanup()
      return True
    
    logging.error('unknown filterMode.')
    return False

      
ARG_OUTPUTTYPE = 'outputType'
ARG_FILTERMODE = 'filterMode'
ARG_OFFSETFIRSTSHRINKWRAP = 'offsetFirstShrinkwrap'
ARG_SPACINGFIRSTREMESH = 'spacingFirstRemesh'
ARG_SPACINGSECONDREMESH = 'spacingSecondRemesh'
ARG_ITERATIONSFIRSTSHRINKWRAP = 'iterationsFirstShrinkwrap'
ARG_ITERATIONSSECONDSHRINKWRAP = 'iterationsSecondShrinkwrap'
ARG_RAYCASTSEARCHEDGELENGTH = 'raycastSearchEdgeLength'
ARG_RAYCASTOUTPUTEDGELENGTH = 'raycastOutputEdgeLength'
ARG_RAYCASTMAXHITDISTANCE = 'raycastMaxHitDistance'
ARG_RAYCASTMAXLENGTH = 'raycastMaxLength'
ARG_RAYCASTMINLENGTH = 'raycastMinLength'
ARG_MAXMODELSDISTANCE = 'maxModelDistance'
ARG_SOLIDIFICATIONTHICKNESS = 'solidificationThickness'
ARG_SMOOTHINGFACTOR = 'smoothingFactor'

OUTPUT_MODEL = 'MODEL'
OUTPUT_SEGMENTATION = 'SEGMENTATION'

MODE_CONVEXHULL = 'CONVEXHULL'
MODE_RAYCASTS = 'RAYCASTS'
MODE_DEEPHULL = 'DEEPHULL'
MODE_NONMANIFOLD = 'NONMANIFOLD'
MODE_SOLIDIFIED = 'SOLIDIFIED'
