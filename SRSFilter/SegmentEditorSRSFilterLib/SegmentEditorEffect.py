import os
import vtk, qt, ctk, slicer
import logging
from SegmentEditorEffects import *
import numpy as np
import math

class SegmentEditorEffect(AbstractScriptedSegmentEditorEffect):
  """This effect uses Watershed algorithm to partition the input volume"""

  def __init__(self, scriptedEffect):
    scriptedEffect.name = 'SRS-Filter'
    scriptedEffect.perSegment = False # this effect operates on all segments at once (not on a single selected segment)
    AbstractScriptedSegmentEditorEffect.__init__(self, scriptedEffect)

    self.segment2DFillOpacity = None
    self.segment2DOutlineOpacity = None
    self.previewedSegmentID = None
    self.filteredImageData = None

    self.timer = qt.QTimer()
    self.previewState = 0
    self.previewStep = 1
    self.previewSteps = 5
    self.timer.connect('timeout()', self.preview)

    self.previewPipelines = {}
    self.setupPreviewDisplay()

    self.logic = SRSFilterLogic(scriptedEffect)
    self.logic.logCallback = self.addLog

    # parameters
    self.parameters = []
    self.offsetFirstShrinkwrapSlider = None
    self.parameters.append({'slider':self.offsetFirstShrinkwrapSlider, 'max':50.0, 'min':1.0,'default':15.0, 'interval':0.1, 'name':'Offset First Shrinkwrap:', 'id':'OffsetFirstShrinkwrap','tooltip':''})
    self.spacingFirstRemeshSlider = None
    self.parameters.append({'slider':self.spacingFirstRemeshSlider, 'max':50.0, 'min':0.1,'default':10.0, 'interval':0.01, 'name':'Spacing First Remesh:', 'id':'SpacingFirstRemesh','tooltip':''})
    self.iterationsFirstShrinkwrapSlider = None
    self.parameters.append({'slider':self.iterationsFirstShrinkwrapSlider, 'max':10, 'min':1,'default':3, 'interval':1, 'name':'Iterations First Shrinkwrap:', 'id':'IterationsFirstShrinkwrap','tooltip':''})
    self.iterationsSecondShrinkwrapSlider = None
    self.parameters.append({'slider':self.iterationsSecondShrinkwrapSlider, 'max':10, 'min':1,'default':5, 'interval':1, 'name':'Iterations Second Shrinkwrap:', 'id':'IterationsSecondShrinkwrap','tooltip':''})
    self.raycastSearchEdgeLengthSlider = None
    self.parameters.append({'slider':self.raycastSearchEdgeLengthSlider, 'max':100.0, 'min':0.1,'default':20.0, 'interval':0.1, 'name':'Raycast Search Edge Length:', 'id':'RaycastSearchEdgeLength','tooltip':''})
    self.raycastOutputEdgeLengthSlider = None
    self.parameters.append({'slider':self.raycastOutputEdgeLengthSlider, 'max':100.0, 'min':0.1,'default':2.0, 'interval':0.1, 'name':'Raycast Output Edge Length:', 'id':'RaycastOutputEdgeLength','tooltip':''})
    self.raycastMaxHitDistanceSlider = None
    self.parameters.append({'slider':self.raycastMaxHitDistanceSlider, 'max':50.0, 'min':0.1,'default':2.0, 'interval':0.01, 'name':'Raycast Max. Hit Distance:', 'id':'RaycastMaxHitDistance','tooltip':''})
    self.raycastMaxLengthSlider = None
    self.parameters.append({'slider':self.raycastMaxLengthSlider, 'max':1000.0, 'min':0.1,'default':100.0, 'interval':0.1, 'name':'Raycast Max. Length:', 'id':'RaycastMaxLength','tooltip':''})
    self.raycastMinLengthSlider = None
    self.parameters.append({'slider':self.raycastMinLengthSlider, 'max':1000.0, 'min':0.0,'default':0.0, 'interval':0.1, 'name':'Raycast Min. Length:', 'id':'RaycastMinLength','tooltip':''})
    self.maxModelsDistanceSlider = None
    self.parameters.append({'slider':self.maxModelsDistanceSlider, 'max':10, 'min':0.01,'default':0.5, 'interval':0.01, 'name':'Max. Models Distance:', 'id':'MaxModelsDistance','tooltip':''})
    self.solidificationThicknessSlider = None
    self.parameters.append({'slider':self.solidificationThicknessSlider, 'max':20.0, 'min':0.1,'default':1.5, 'interval':0.1, 'name':'Solidification Thickness:', 'id':'SolidificationThickness','tooltip':''})
    
    # filter modes
    self.filterModeTypeMap = {}
    self.filterModes = []
    self.surfaceButton = None
    self.filterModes.append({'button':self.surfaceButton, 'name':'Surface', 'id':'SURFACE', 'default':True})
    self.hullShallowButton = None
    self.filterModes.append({'button':self.hullShallowButton, 'name':'Hull Shallow', 'id':'SHALLOW','default':False})
    self.raycastResultButton = None
    self.filterModes.append({'button':self.raycastResultButton, 'name':'Raycast Result', 'id':'RAYCAST','default':False})
    self.hullDeepButton = None
    self.filterModes.append({'button':self.hullDeepButton, 'name':'Hull Deep', 'id':'DEEP','default':False})
    self.nonmanifoldModelButton = None
    self.filterModes.append({'button':self.nonmanifoldModelButton, 'name':'Solid Model', 'id':'NONSOLID','default':False})
    self.manifoldModelButton = None
    self.filterModes.append({'button':self.manifoldModelButton, 'name':'Non-Manifold Model', 'id':'SOLID','default':False})

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
    return """bla"""

  def activate(self):
    pass

  def deactivate(self):
    self.restorePreviewedSegmentTransparency()

    # Clear preview pipeline and stop timer
    self.clearPreviewDisplay()
    self.timer.stop()
    self.filteredOrientedImageData = None

  def setCurrentSegmentTransparent(self):
    """Save current segment opacity and set it to zero
    to temporarily hide the segment so that threshold preview
    can be seen better.
    It also restores opacity of previously previewed segment.
    Call restorePreviewedSegmentTransparency() to restore original
    opacity.
    """
    segmentationNode = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    if not segmentationNode:
      return
    displayNode = segmentationNode.GetDisplayNode()
    if not displayNode:
      return
    segmentID = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()

    if segmentID == self.previewedSegmentID:
      # already previewing the current segment
      return

    # If an other segment was previewed before, restore that.
    if self.previewedSegmentID:
      self.restorePreviewedSegmentTransparency()

    # Make current segment fully transparent
    if segmentID:
      self.segment2DFillOpacity = displayNode.GetSegmentOpacity2DFill(segmentID)
      self.segment2DOutlineOpacity = displayNode.GetSegmentOpacity2DOutline(segmentID)
      self.previewedSegmentID = segmentID
      displayNode.SetSegmentOpacity2DFill(segmentID, 0)
      displayNode.SetSegmentOpacity2DOutline(segmentID, 0)
  
  def restorePreviewedSegmentTransparency(self):
    """Restore previewed segment's opacity that was temporarily
    made transparen by calling setCurrentSegmentTransparent()."""
    segmentationNode = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    if not segmentationNode:
      return
    displayNode = segmentationNode.GetDisplayNode()
    if not displayNode:
      return
    if not self.previewedSegmentID:
      # already previewing the current segment
      return
    displayNode.SetSegmentOpacity2DFill(self.previewedSegmentID, self.segment2DFillOpacity)
    displayNode.SetSegmentOpacity2DOutline(self.previewedSegmentID, self.segment2DOutlineOpacity)
    self.previewedSegmentID = None

  def setupOptionsFrame(self):
    
    self.updatePreviewButton = qt.QPushButton("Update Preview")
    self.updatePreviewButton.objectName = self.__class__.__name__ + 'Update Preview'
    self.updatePreviewButton.setToolTip("Fill selected segment in regions that are in the specified intensity range.")
    self.scriptedEffect.addOptionsWidget(self.updatePreviewButton)

    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.objectName = self.__class__.__name__ + 'Apply'
    self.applyButton.setToolTip("Fill selected segment in regions that are in the specified intensity range.")
    self.scriptedEffect.addOptionsWidget(self.applyButton)

    self.applyButton.connect('clicked()', self.onApply)
    self.updatePreviewButton.connect('clicked()', self.onUpdatePreview)

    # setup parameters
    for param in self.parameters:
      param['slider'] = slicer.qMRMLSliderWidget()
      param['slider'].setMRMLScene(slicer.mrmlScene)
      #param['slider'].quantity = "length" # get unit, precision, etc. from MRML unit node
      param['slider'].minimum = param['min']
      param['slider'].maximum = param['max']
      param['slider'].tickInterval = param['interval']
      param['slider'].value = param['default']
      param['slider'].setToolTip(param['tooltip'])
      self.scriptedEffect.addLabeledOptionsWidget(param['name'], param['slider'])
      param['slider'].connect('valueChanged(double)', self.updateMRMLFromGUI)

    # setup filter modes
    filterModeLayout = qt.QGridLayout()
    columnCount = 0
    rowCount = 0
    self.filterModeButtons = []
    for mode in self.filterModes:
      mode['button'] = qt.QRadioButton(mode['name'])
      self.filterModeButtons.append(mode['button'])
      self.filterModeTypeMap[mode['button']] = mode['id']
      mode['button'].connect('toggled(bool)', self.updateMRMLFromGUI)
      
      filterModeLayout.addWidget(mode['button'], rowCount, columnCount)
      
      if columnCount >= 2:
        columnCount = 0
        rowCount += 1
      else:
        columnCount += 1

    self.scriptedEffect.addLabeledOptionsWidget("Filter Mode:", filterModeLayout)


  def createCursor(self, widget):
    # Turn off effect-specific cursor for this effect
    return slicer.util.mainWindow().cursor

  def layoutChanged(self):
    self.setupPreviewDisplay()

  def processInteractionEvents(self, callerInteractor, eventId, viewWidget):
    return False # For the sake of example

  def processViewNodeEvents(self, callerViewNode, eventId, viewWidget):
    pass # For the sake of example

  def setMRMLDefaults(self):
    for param in self.parameters:
      self.scriptedEffect.setParameterDefault(param['id'], param['default'])

    self.scriptedEffect.setParameterDefault("Filtermode", next(item for item in self.filterModes if item['default'] == True)['id'])

  def updateGUIFromMRML(self):
    for param in self.parameters:
      value = self.scriptedEffect.doubleParameter(param['id'])
      wasBlocked = param['slider'].blockSignals(True)
      param['slider'].value = abs(value)
      param['slider'].blockSignals(wasBlocked)
    
    filterModeName = self.scriptedEffect.parameter("Filtermode")
    filterModeButton = list(self.filterModeTypeMap.keys())[list(self.filterModeTypeMap.values()).index(filterModeName)]
    filterModeButton.setChecked(True)

  def updateMRMLFromGUI(self):
    for param in self.parameters:
      self.scriptedEffect.setParameter(param['id'], param['slider'].value)

    for button in self.filterModeTypeMap:
      if button.isChecked():
        self.scriptedEffect.setParameter("Filtermode", self.filterModeTypeMap[button])
  

  #
  # Effect specific methods (the above ones are the API methods to override)
  #

  def onUpdatePreview(self):

    # run SRS-Filter
    self.createFilteredImageData()

    # Setup and start preview pulse
    self.setCurrentSegmentTransparent()
    self.setupPreviewDisplay()
    self.timer.start(200)

  def onApply(self):

    self.scriptedEffect.saveStateForUndo()

    # Apply changes
    # check if self.filteredOrientedImageData is empty

    self.scriptedEffect.modifySelectedSegmentByLabelmap(self.filteredOrientedImageData, slicer.qSlicerSegmentEditorAbstractEffect.ModificationModeSet)
    self.filteredOrientedImageData = None

    # De-select effect
    self.scriptedEffect.selectEffect("")

  def createFilteredImageData(self):
    
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

    # Get master volume image data
    import vtkSegmentationCorePython

    # Get modifier labelmap
    selectedSegmentLabelmap = self.scriptedEffect.selectedSegmentLabelmap()

    #TODO: check if segment is not empty

    self.logic.srsFilter(self.scriptedEffect.defaultModifierLabelmap(), self.filteredOrientedImageData)

    slicer.util.showStatusMessage('')
    qt.QApplication.restoreOverrideCursor()
    
    #self.filteredOrientedImageData.DeepCopy(thresh.GetOutput())
    

  def clearPreviewDisplay(self):
    for sliceWidget, pipeline in self.previewPipelines.items():
      self.scriptedEffect.removeActor2D(sliceWidget, pipeline.actor)
    self.previewPipelines = {}

  def setupPreviewDisplay(self):
    # Clear previous pipelines before setting up the new ones
    self.clearPreviewDisplay()

    layoutManager = slicer.app.layoutManager()
    if layoutManager is None:
      return

    # Add a pipeline for each 2D slice view
    for sliceViewName in layoutManager.sliceViewNames():
      sliceWidget = layoutManager.sliceWidget(sliceViewName)
      if not self.scriptedEffect.segmentationDisplayableInView(sliceWidget.mrmlSliceNode()):
        continue
      renderer = self.scriptedEffect.renderer(sliceWidget)
      if renderer is None:
        logging.error("setupPreviewDisplay: Failed to get renderer!")
        continue

      # Create pipeline
      pipeline = PreviewPipeline()
      self.previewPipelines[sliceWidget] = pipeline

      # Add actor
      self.scriptedEffect.addActor2D(sliceWidget, pipeline.actor)

  def preview(self):
    """Set values to pipeline"""

    import vtkSegmentationCorePython

    opacity = 0.5 + self.previewState / (2. * self.previewSteps)
    
    if self.filteredOrientedImageData is None:
      logging.error('SRS-Filter not applied yet.')
      return

    # Get color of edited segment
    segmentationNode = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    if not segmentationNode:
      # scene was closed while preview was active
      return
    displayNode = segmentationNode.GetDisplayNode()
    if displayNode is None:
      logging.error("preview: Invalid segmentation display node!")
      color = [0.5,0.5,0.5]
    segmentID = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()

    # Make sure we keep the currently selected segment hidden (the user may have changed selection)
    if segmentID != self.previewedSegmentID:
      self.setCurrentSegmentTransparent()

    r,g,b = segmentationNode.GetSegmentation().GetSegment(segmentID).GetColor()

    img = vtkSegmentationCorePython.vtkOrientedImageData()
    img.DeepCopy(self.filteredOrientedImageData)
    img.SetImageToWorldMatrix(vtk.vtkMatrix4x4())

    # Set values to pipelines
    for sliceWidget in self.previewPipelines:
      pipeline = self.previewPipelines[sliceWidget]
      pipeline.lookupTable.SetTableValue(1,  r, g, b,  opacity)
      sliceLogic = sliceWidget.sliceLogic()
      backgroundLogic = sliceLogic.GetBackgroundLayer()
      
      #pipeline.thresholdFilter.ThresholdBetween(min, max)
      ##slice here
      #sliceNode = sliceWidget.mrmlSliceNode()
      reslice = vtk.vtkImageReslice()
      reslice.SetInputData(img)
      reslice.SetOutputDimensionality(2)
      reslice.SetInterpolationModeToLinear()
      reslice.SetResliceTransform(backgroundLogic.GetReslice().GetResliceTransform())
      
      pipeline.colorMapper.SetInputConnection(reslice.GetOutputPort())
      pipeline.actor.VisibilityOn()
      sliceWidget.sliceView().scheduleRender()

    self.previewState += self.previewStep
    if self.previewState >= self.previewSteps:
      self.previewStep = -1
    if self.previewState <= 0:
      self.previewStep = 1


  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update



class SRSFilterLogic(object):

  def __init__(self, scriptedEffect):
    self.scriptedEffect = scriptedEffect
    self.logCallback = None

  def srsFilter(self, inputImageData, outputImageData):

    def polydataToImagedata(polydata):

      outputImageData.DeepCopy(inputImageData)
      
      pol2stenc = vtk.vtkPolyDataToImageStencil()
      pol2stenc.SetInputData(polydata)
      pol2stenc.SetOutputOrigin(inputImageData.GetOrigin())
      pol2stenc.SetOutputSpacing(inputImageData.GetSpacing())
      pol2stenc.SetOutputWholeExtent(inputImageData.GetExtent())
      pol2stenc.Update()

      imgstenc = vtk.vtkImageStencil()
      imgstenc.SetInputData(outputImageData)
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

      outputImageData.DeepCopy(revimgstenc.GetOutput())
      
      return outputImageData

    def smoothPolydata(polydata):
      decimator = vtk.vtkDecimatePro()
      decimator.SetInputData(polydata)
      decimator.SetFeatureAngle(60)
      decimator.SplittingOff()
      decimator.PreserveTopologyOn()
      decimator.SetMaximumError(1)
      decimator.SetTargetReduction(0.25)
      decimator.ReleaseDataFlagOff()
      decimator.Update()

      smootherSinc = vtk.vtkWindowedSincPolyDataFilter()
      smootherSinc.SetPassBand(0.1)
      smootherSinc.SetInputConnection(decimator.GetOutputPort())
      smootherSinc.SetNumberOfIterations(10)
      smootherSinc.FeatureEdgeSmoothingOff()
      smootherSinc.BoundarySmoothingOff()
      smootherSinc.ReleaseDataFlagOn()
      smootherSinc.Update()

      return smootherSinc.GetOutput()

    OFFSETFIRSTSHRINKWRAP = self.scriptedEffect.doubleParameter('OffsetFirstShrinkwrap')
    SPACINGFIRSTREMESH = self.scriptedEffect.doubleParameter('SpacingFirstRemesh')
    ITERATIONSFIRSTSHRINKWRAP = int(self.scriptedEffect.doubleParameter('IterationsFirstShrinkwrap'))
    ITERATIONSSECONDSHRINKWRAP = int(self.scriptedEffect.doubleParameter('IterationsSecondShrinkwrap'))
    RAYCASTSEARCHEDGELENGTH = self.scriptedEffect.doubleParameter('RaycastSearchEdgeLength')
    RAYCASTOUTPUTEDGELENGTH = self.scriptedEffect.doubleParameter('RaycastOutputEdgeLength')
    RAYCASTMAXHITDISTANCE = self.scriptedEffect.doubleParameter('RaycastMaxHitDistance')
    RAYCASTMAXLENGTH = self.scriptedEffect.doubleParameter('RaycastMaxLength')
    RAYCASTMINLENGTH = self.scriptedEffect.doubleParameter('RaycastMinLength')
    MAXMODELSDISTANCE = self.scriptedEffect.doubleParameter('MaxModelsDistance')
    THICKNESS = self.scriptedEffect.doubleParameter('SolidificationThickness')
    FILTERMODE = self.scriptedEffect.parameter('Filtermode')

    # OFFSETFIRSTSHRINKWRAP = 15# self.scriptedEffect.doubleParameter('OffsetFirstShrinkwrap')
    # SPACINGFIRSTREMESH = 10#self.scriptedEffect.doubleParameter('SpacingFirstRemesh')
    # ITERATIONSFIRSTSHRINKWRAP = 3#int(self.scriptedEffect.doubleParameter('IterationsFirstShrinkwrap'))
    # ITERATIONSSECONDSHRINKWRAP = 5#int(self.scriptedEffect.doubleParameter('IterationsSecondShrinkwrap'))
    # RAYCASTSEARCHEDGELENGTH = 20#self.scriptedEffect.doubleParameter('RaycastSearchEdgeLength')
    # RAYCASTOUTPUTEDGELENGTH = 2#self.scriptedEffect.doubleParameter('RaycastOutputEdgeLength')
    # RAYCASTMAXHITDISTANCE = 2#self.scriptedEffect.doubleParameter('RaycastMaxHitDistance')
    # RAYCASTMAXLENGTH = 100#self.scriptedEffect.doubleParameter('RaycastMaxLength')
    # RAYCASTMINLENGTH = 0#self.scriptedEffect.doubleParameter('RaycastMinLength')
    # MAXMODELSDISTANCE = 0.5#elf.scriptedEffect.doubleParameter('MaxModelsDistance')
    # THICKNESS = 1.5#self.scriptedEffect.doubleParameter('SolidificationThickness')
    # FILTERMODE = 'SHALLOW'#self.scriptedEffect.parameter('Filtermode')


    # create model from segmentation
    threshold = vtk.vtkImageThreshold()
    threshold.SetInputData(inputImageData)
    threshold.ThresholdBetween(0,0)
    threshold.ReplaceOutOn()
    threshold.ReplaceInOn()
    threshold.SetInValue(0)
    threshold.SetOutValue(1)
    threshold.Update()

    inputDiscreteCubes = vtk.vtkDiscreteMarchingCubes()
    inputDiscreteCubes.SetInputData(threshold.GetOutput())
    inputDiscreteCubes.GenerateValues(1,0,0)
    inputDiscreteCubes.Update()

    #region create sphere
    bounds = np.array([0]*6)
    inputDiscreteCubes.GetOutput().GetBounds(bounds)
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
    cellLocator.SetDataSet(inputDiscreteCubes.GetOutput())
    cellLocator.BuildLocator()

    for x in range(ITERATIONSFIRSTSHRINKWRAP):
      
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

        if OFFSETFIRSTSHRINKWRAP > 0 and vectorLength > 0.01:
          newLocation = closestPoint - ((vector/vectorLength) * OFFSETFIRSTSHRINKWRAP)
        else:
          newLocation = closestPoint
        
        points.SetPoint(i, newLocation)
      
      shrinkModelPD.SetPoints(points)

      if x == (ITERATIONSFIRSTSHRINKWRAP - 1):
        if self.logCallback: self.logCallback('First shrinkwrap completed.')
        break

      # remesh
      
      whiteImage = vtk.vtkImageData()
      bounds = [0]*6
      shrinkModelPD.GetBounds(bounds)

      spacing = [SPACINGFIRSTREMESH]*3
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
      pol2stenc.SetInputData(shrinkModelPD)

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

      shrinkModelPD.DeepCopy(reverse.GetOutput())
      if self.logCallback: self.logCallback('First shrinkwrap: %s/%s.' %(x+1, ITERATIONSFIRSTSHRINKWRAP))
    
    if FILTERMODE == 'SHALLOW':
      return polydataToImagedata(shrinkModelPD)

    #endregion

    #region Raycast

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

      if max(edgeLength) > RAYCASTSEARCHEDGELENGTH:
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
    adapt.SetMaximumEdgeLength(RAYCASTOUTPUTEDGELENGTH)
    adapt.SetMaximumTriangleArea(vtk.VTK_INT_MAX)
    adapt.SetMaximumNumberOfPasses(vtk.VTK_INT_MAX)
    adapt.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(adapt.GetOutput())
    clean.Update()

    shrinkModelPD.DeepCopy(clean.GetOutput())

    if largeCellIds.GetNumberOfIds() > 0 and RAYCASTMAXLENGTH > 0.0:

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
              vector = normal * RAYCASTMAXLENGTH # max Length greater 0, checked above

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

              loc_new = np.array(glo) - (normal * 0.5)
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
          if vert_location_dict[i][2] > RAYCASTMINLENGTH:
            
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
                  if distance < RAYCASTMAXHITDISTANCE:
                    shrinkModelPD.GetPoints().SetPoint(i, vert_location_dict[i][1])
                    pointChanged = True
    if self.logCallback: self.logCallback('Raycast completed.')

    if FILTERMODE == 'RAYCAST':
      return polydataToImagedata(shrinkModelPD)
    
    #endregion

    #region Shrinkwrap

    for x in range(ITERATIONSSECONDSHRINKWRAP):
      
      # remesh
      
      whiteImage = vtk.vtkImageData()
      bounds = [0]*6
      shrinkModelPD.GetBounds(bounds)

      #spacing = [SPACINGSECONDREMESH]*3
      spacing = inputImageData.GetSpacing()
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
      pol2stenc.SetInputData(shrinkModelPD)

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

      shrinkModelPD.DeepCopy(reverse.GetOutput())

      if x == (ITERATIONSSECONDSHRINKWRAP - 1):
        if self.logCallback: self.logCallback('Second shrinkwrap completed.')
        break

      # shrinkwrap

      smoothFilter = vtk.vtkSmoothPolyDataFilter()
      smoothFilter.SetInputData(0, shrinkModelPD)
      smoothFilter.SetInputData(1, inputDiscreteCubes.GetOutput())
      smoothFilter.Update()
      
      shrinkModelPD.DeepCopy(smoothFilter.GetOutput())

      if self.logCallback: self.logCallback('Second shrinkwrap: %s/%s.' %(x+1, ITERATIONSSECONDSHRINKWRAP))
    
    if FILTERMODE == 'DEEP':
      return polydataToImagedata(shrinkModelPD)

    #endregion

    # region Remove Caps

    # implicit distance, add point ids with larger distance to ids
    implicitDistance = vtk.vtkImplicitPolyDataDistance()
    implicitDistance.SetInput(inputDiscreteCubes.GetOutput())
    
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

        if abs(distance) > MAXMODELSDISTANCE:
          nonsolidPolyData.DeleteCell(c)
          break

    nonsolidPolyData.RemoveDeletedCells()
    shrinkModelPD.DeepCopy(nonsolidPolyData)
    if self.logCallback: self.logCallback('Caps removed.')

    if FILTERMODE == 'NONSOLID':
      
      modelsLogic = slicer.modules.models.logic()
      modelNode = modelsLogic.AddModel(smoothPolydata(shrinkModelPD))
      seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode().GetSegmentation().GetSegment(self.scriptedEffect.parameterSetNode().GetSelectedSegmentID())
      modelNode.GetDisplayNode().SetColor(seg.GetColor())
      modelNode.SetName(seg.GetName())
      outputImageData.DeepCopy(inputImageData)
      
      return outputImageData

    #endregion

    #region Solidification

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

    polydata = vtk.vtkPolyData()
    polydata.DeepCopy(normals.GetOutput())
    numberOfPoints = polydata.GetNumberOfPoints()

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
      polydata.GetPointCells(pointID, cellIDs)
      normalsArray = []

      
      # ilterate through all cells/faces which contain point
      for id in xrange(cellIDs.GetNumberOfIds()):
        # faceData = []
        n = []
        n.append(polydata.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3))
        n.append(polydata.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3 + 1))
        n.append(polydata.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(id)*3 + 2))

        normalsArray.append(np.array(n) * (-1))

      # calculate position of new vert
      dir_vec = np.zeros(3)
      
      for n in normalsArray:
        dir_vec = dir_vec + np.array(n)

      dir_vec_norm = dir_vec / np.linalg.norm(dir_vec)
      proj_length = np.dot(dir_vec_norm, np.array(normalsArray[0]))
      dir_vec_finallenght = dir_vec_norm * proj_length
      vertex_neu = np.array(polydata.GetPoint(pointID)) + (dir_vec_finallenght * THICKNESS)
      
      # append point
      addingPoints.append(vertex_neu)

    for cellID in range(polydata.GetNumberOfCells()):
      pointIDs = vtk.vtkIdList()
      polydata.GetCellPoints(cellID, pointIDs)

      newPointIDs = vtk.vtkIdList()
      for i in reversed(range(pointIDs.GetNumberOfIds())):
        newPointIDs.InsertNextId(int(pointIDs.GetId(i) + numberOfPoints))

      addingPolys.append(newPointIDs)

    doubleSurfacePoints = vtk.vtkPoints()
    doubleSurfacePolys = vtk.vtkCellArray()

    doubleSurfacePoints.DeepCopy(polydata.GetPoints())
    doubleSurfacePolys.DeepCopy(polydata.GetPolys())

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

    shrinkModelPD.DeepCopy(triangleFilter.GetOutput())
    if self.logCallback: self.logCallback('Model solidified.')

    if FILTERMODE == 'SOLID':
      modelsLogic = slicer.modules.models.logic()
      modelNode = modelsLogic.AddModel(smoothPolydata(shrinkModelPD))
      seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode().GetSegmentation().GetSegment(self.scriptedEffect.parameterSetNode().GetSelectedSegmentID())
      modelNode.GetDisplayNode().SetColor(seg.GetColor())
      modelNode.SetName(seg.GetName())
      outputImageData.DeepCopy(inputImageData)
      
      return outputImageData

    return polydataToImagedata(shrinkModelPD)
    #endregion
      

#
# PreviewPipeline
#
class PreviewPipeline(object):
  """ Visualization objects and pipeline for each slice view for threshold preview
  """

  def __init__(self):
    self.lookupTable = vtk.vtkLookupTable()
    self.lookupTable.SetRampToLinear()
    self.lookupTable.SetNumberOfTableValues(2)
    self.lookupTable.SetTableRange(0, 1)
    self.lookupTable.SetTableValue(0,  0, 0, 0,  0)
    self.colorMapper = vtk.vtkImageMapToRGBA()
    self.colorMapper.SetOutputFormatToRGBA()
    self.colorMapper.SetLookupTable(self.lookupTable)
    self.thresholdFilter = vtk.vtkImageThreshold()
    self.thresholdFilter.SetInValue(1)
    self.thresholdFilter.SetOutValue(0)
    self.thresholdFilter.SetOutputScalarTypeToUnsignedChar()

    # Feedback actor
    self.mapper = vtk.vtkImageMapper()
    self.dummyImage = vtk.vtkImageData()
    self.dummyImage.AllocateScalars(vtk.VTK_UNSIGNED_INT, 1)
    self.mapper.SetInputData(self.dummyImage)
    self.actor = vtk.vtkActor2D()
    self.actor.VisibilityOff()
    self.actor.SetMapper(self.mapper)
    self.mapper.SetColorWindow(255)
    self.mapper.SetColorLevel(128)

    # Setup pipeline
    #self.colorMapper.SetInputConnection(self.thresholdFilter.GetOutputPort())
    self.mapper.SetInputConnection(self.colorMapper.GetOutputPort())

