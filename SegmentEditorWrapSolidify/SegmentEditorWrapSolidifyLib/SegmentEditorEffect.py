import os
import vtk, qt, ctk, slicer
import logging
from SegmentEditorEffects import *
import numpy as np
import math
import vtkSegmentationCorePython

class SegmentEditorEffect(AbstractScriptedSegmentEditorEffect):
  """This effect uses shrinkwrap, raycasting, remesh, and solidifying algorithms to filter the surface from the input segmentation"""

  def __init__(self, scriptedEffect):
    scriptedEffect.name = 'Wrap Solidify'
    scriptedEffect.perSegment = True # this effect operates on all segments at once (not on a single selected segment)
    AbstractScriptedSegmentEditorEffect.__init__(self, scriptedEffect)

    self.logic = WrapSolidifyLogic(scriptedEffect)
    self.logic.logCallback = self.addLog


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
    # TODO: adapt
    return """<html>This Filter results in an even surface of the current segment. It is using a combination of shrinkwrapping, projection and solidification algorithms.<br>
    For further information, license, disclaimers and possible research partnerships visit <a href="https://github.com/sebastianandress/Slicer-SurfaceWrapSolidify">this</a> github repository.
    </html>"""

  def activate(self):
    pass

  def deactivate(self):
    self.cleanup()
  
  def cleanup(self):
    pass


  def setupOptionsFrame(self):

    # smoothing factor

    self.smoothingFactorSlider = slicer.qMRMLSliderWidget()
    self.smoothingFactorSlider.setMRMLScene(slicer.mrmlScene)
    self.smoothingFactorSlider.minimum = 0
    self.smoothingFactorSlider.maximum = 1
    self.smoothingFactorSlider.singleStep = 0.1
    self.smoothingFactorSlider.value = DEFAULT_SMOOTHINGFACTOR
    self.smoothingFactorSlider.setToolTip('Smoothing Factor used for operation, as the a surface representation of the segmentation will be used for this algorithm.')
    self.smoothingFactorSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.scriptedEffect.addLabeledOptionsWidget('Smoothing Factor: ', self.smoothingFactorSlider)


    self.scriptedEffect.addOptionsWidget(qt.QLabel(''))

    # carve out cavities

    self.cavitiesCheckBox = qt.QCheckBox('Carve out Cavities')
    self.cavitiesCheckBox.setChecked(DEFAULT_CARVECAVITIES)
    self.cavitiesCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)
    self.scriptedEffect.addOptionsWidget(self.cavitiesCheckBox)

    self.cavitiesDiameterSlider = slicer.qMRMLSliderWidget()
    self.cavitiesDiameterSlider.setMRMLScene(slicer.mrmlScene)
    self.cavitiesDiameterSlider.setEnabled(False)
    self.cavitiesDiameterSlider.minimum = 0.1
    self.cavitiesDiameterSlider.maximum = 100.0
    self.cavitiesDiameterSlider.singleStep = 0.1
    self.cavitiesDiameterSlider.value = DEFAULT_CAVITIESDIAMETER
    self.cavitiesDiameterSlider.suffix = 'mm'
    self.cavitiesDiameterSlider.setToolTip('Cavities with larger diameter will not be solidified.')
    self.cavitiesDiameterSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.scriptedEffect.addLabeledOptionsWidget('   Minimal Cavities Diameter: ', self.cavitiesDiameterSlider)

    self.cavitiesDepthSlider = slicer.qMRMLSliderWidget()
    self.cavitiesDepthSlider.setMRMLScene(slicer.mrmlScene)
    self.cavitiesDepthSlider.setEnabled(False)
    self.cavitiesDepthSlider.minimum = 0.1
    self.cavitiesDepthSlider.maximum = 1000.0
    self.cavitiesDepthSlider.singleStep = 0.1
    self.cavitiesDepthSlider.value = DEFAULT_CAVITIESDEPTH
    self.cavitiesDepthSlider.suffix = 'mm'
    self.cavitiesDepthSlider.setToolTip('Deeper cavities will not be solidified.')
    self.cavitiesDepthSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.scriptedEffect.addLabeledOptionsWidget('   Minimal Cavities Depth: ', self.cavitiesDepthSlider)

    # create shell

    self.createShellCheckBox = qt.QCheckBox('Create Shell')
    self.createShellCheckBox.setChecked(DEFAULT_CREATESHELL)
    self.createShellCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)
    self.scriptedEffect.addOptionsWidget(self.createShellCheckBox)

    self.shellThicknessSlider = slicer.qMRMLSliderWidget()
    self.shellThicknessSlider.setMRMLScene(slicer.mrmlScene)
    self.shellThicknessSlider.setEnabled(False)
    self.shellThicknessSlider.minimum = -0.1
    self.shellThicknessSlider.maximum = 20.0
    self.shellThicknessSlider.singleStep = 0.1
    self.shellThicknessSlider.value = DEFAULT_SHELLTHICKNESS
    self.shellThicknessSlider.suffix = 'mm'
    self.shellThicknessSlider.setToolTip('Thickness of the output shell.\nCAVE: If this smaller than the spacing of the input segmentation, it might appear punctured in the output. Please select "Model" as output type then. If <0, a nonmanifold mesh gets created, Segmentation will fail.')
    self.shellThicknessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.scriptedEffect.addLabeledOptionsWidget('   Output Shell Thickness: ', self.shellThicknessSlider)


    self.scriptedEffect.addOptionsWidget(qt.QLabel(''))

    # output types

    self.outputTypeLayout = qt.QVBoxLayout()
    self.scriptedEffect.addLabeledOptionsWidget('Output Type: ', self.outputTypeLayout)
    self.outputTypeGroup = qt.QButtonGroup()

    self.outputSegmentationRadioButton = qt.QRadioButton('Segmentation')
    self.outputTypeGroup.addButton(self.outputSegmentationRadioButton)
    self.outputTypeLayout.addWidget(self.outputSegmentationRadioButton)
    if DEFAULT_OUTPUTTYPE == OUTPUT_SEGMENTATION:
      self.outputSegmentationRadioButton.setChecked(True)


    modelFrame = qt.QFrame()
    modelLayout = qt.QHBoxLayout()
    modelLayout.setContentsMargins(0,0,0,0)
    modelFrame.setLayout(modelLayout)
    self.outputTypeLayout.addWidget(modelFrame)

    self.outputModelRadioButton = qt.QRadioButton('Model: ')
    self.outputTypeGroup.addButton(self.outputModelRadioButton)
    modelFrame.layout().addWidget(self.outputModelRadioButton)
    if DEFAULT_OUTPUTTYPE == OUTPUT_MODEL:
      self.outputModelRadioButton.setChecked(True)

    self.outputTypeGroup.connect('buttonClicked(int)', self.updateMRMLFromGUI)

    self.outputModelSelector = slicer.qMRMLNodeComboBox()
    self.outputModelSelector.nodeTypes = ["vtkMRMLModelNode"]
    self.outputModelSelector.selectNodeUponCreation = False
    self.outputModelSelector.addEnabled = False
    self.outputModelSelector.renameEnabled = True
    self.outputModelSelector.removeEnabled = True
    self.outputModelSelector.noneEnabled = True
    self.outputModelSelector.noneDisplay = 'Create new Model'
    self.outputModelSelector.showHidden = False
    self.outputModelSelector.showChildNodeTypes = False
    self.outputModelSelector.baseName = "Model"
    self.outputModelSelector.selectNodeUponCreation = True
    self.outputModelSelector.setMRMLScene(slicer.mrmlScene)
    modelFrame.layout().addWidget(self.outputModelSelector)
    self.outputModelSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.updateMRMLFromGUI)


    # advanced options

    self.scriptedEffect.addOptionsWidget(qt.QLabel(''))

    advancedSettingsFrame = ctk.ctkCollapsibleGroupBox()
    advancedSettingsFrame.objectName = 'advancedSettingsFrame'
    advancedSettingsFrame.title = "Advanced Options"
    advancedSettingsFrame.collapsed = True
    advancedSettingsFrame.setLayout(qt.QFormLayout())
    self.scriptedEffect.addOptionsWidget(advancedSettingsFrame)

    # nr of iterations (always last), default 5
    self.iterationsSlider = slicer.qMRMLSliderWidget()
    self.iterationsSlider.setMRMLScene(slicer.mrmlScene)
    self.iterationsSlider.minimum = 1
    self.iterationsSlider.maximum = 20
    self.iterationsSlider.singleStep = 1
    self.iterationsSlider.value = DEFAULT_ITERATIONS
    self.iterationsSlider.suffix = ''
    self.iterationsSlider.setToolTip('Increase this value if output seems not completely converged to input.\nCAVE: Increases computation time.')
    advancedSettingsFrame.layout().addRow('Number of Shrinkwrap Iterations: ', self.iterationsSlider)
    self.iterationsSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    
    # spacing (> 1 or volume min of all spacings)
    self.spacingSlider = slicer.qMRMLSliderWidget()
    self.spacingSlider.setMRMLScene(slicer.mrmlScene)
    self.spacingSlider.minimum = 0.1
    self.spacingSlider.maximum = 10
    self.spacingSlider.singleStep = 0.1
    self.spacingSlider.value = DEFAULT_SPACING
    #TODO: set default depending on input volume
    self.spacingSlider.suffix = 'mm^3'
    self.spacingSlider.setToolTip('Increase this value if output seems terraced.\nCAVE: Drastically increases computation time.')
    advancedSettingsFrame.layout().addRow('Remesh Spacing: ', self.spacingSlider)
    self.spacingSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)

    # shell to input distance
    self.shellDistanceSlider = slicer.qMRMLSliderWidget()
    self.shellDistanceSlider.setMRMLScene(slicer.mrmlScene)
    self.shellDistanceSlider.minimum = -0.1
    self.shellDistanceSlider.maximum = 10
    self.shellDistanceSlider.singleStep = 0.1
    self.shellDistanceSlider.value = DEFAULT_SHELLDISTANCE
    self.shellDistanceSlider.suffix = 'mm'
    self.shellDistanceSlider.setToolTip('Increase this value if output seems punctated. If <0 is selected, no parts will be removed.\nCAVE: Might bridge areas (and therefore for example hide fracture gaps).')
    advancedSettingsFrame.layout().addRow('Shell to Input Distance: ', self.shellDistanceSlider)
    self.shellDistanceSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    
    # Apply Button

    self.scriptedEffect.addOptionsWidget(qt.QLabel(''))
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
    self.scriptedEffect.setParameterDefault(ARG_CARVECAVITIES, DEFAULT_CARVECAVITIES)
    self.scriptedEffect.setParameterDefault(ARG_CAVITIESDIAMETER, DEFAULT_CAVITIESDIAMETER)
    self.scriptedEffect.setParameterDefault(ARG_CAVITIESDEPTH, DEFAULT_CAVITIESDEPTH)

    self.scriptedEffect.setParameterDefault(ARG_CREATESHELL, DEFAULT_CREATESHELL)
    self.scriptedEffect.setParameterDefault(ARG_SHELLTHICKNESS, DEFAULT_SHELLTHICKNESS)

    self.scriptedEffect.setParameterDefault(ARG_OUTPUTTYPE, DEFAULT_OUTPUTTYPE)
    self.scriptedEffect.setParameterDefault(ARG_SMOOTHINGFACTOR, DEFAULT_SMOOTHINGFACTOR)
    self.scriptedEffect.setParameterDefault(ARG_OUTPUTMODELNODE, None)
    
    self.scriptedEffect.setParameterDefault(ARG_ITERATIONS, DEFAULT_ITERATIONS)
    self.scriptedEffect.setParameterDefault(ARG_SPACING, DEFAULT_SPACING)
    self.scriptedEffect.setParameterDefault(ARG_SHELLDISTANCE, DEFAULT_SHELLDISTANCE)

  def updateGUIFromMRML(self):
    for pId, pElement in [
        (ARG_SMOOTHINGFACTOR, self.smoothingFactorSlider),
        (ARG_CAVITIESDIAMETER, self.cavitiesDiameterSlider),
        (ARG_CAVITIESDEPTH, self.cavitiesDepthSlider),
        (ARG_SHELLTHICKNESS, self.shellThicknessSlider),
        (ARG_ITERATIONS, self.iterationsSlider),
        (ARG_SPACING, self.spacingSlider),
        (ARG_SHELLDISTANCE, self.shellDistanceSlider)
      ]:
      value = self.scriptedEffect.doubleParameter(pId)
      wasBlocked = pElement.blockSignals(True)
      pElement.value = value
      pElement.blockSignals(wasBlocked)

    self.cavitiesCheckBox.setChecked(self.scriptedEffect.parameter(ARG_CARVECAVITIES)=='True')
    self.createShellCheckBox.setChecked(self.scriptedEffect.parameter(ARG_CREATESHELL)=='True')
    
    outputTypeID = self.scriptedEffect.parameter(ARG_OUTPUTTYPE)
    if OUTPUT_SEGMENTATION == outputTypeID:
      self.outputSegmentationRadioButton.setChecked(True)
      self.outputModelRadioButton.setChecked(False)
    elif OUTPUT_MODEL == outputTypeID:
      self.outputModelRadioButton.setChecked(True)
      self.outputSegmentationRadioButton.setChecked(False)
    
    self.outputModelSelector.setCurrentNode(slicer.mrmlScene.GetNodeByID(self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE)))

    self.disableOptions()
    self.cleanup()

  def updateMRMLFromGUI(self):
    for pId, pElement in [
        (ARG_SMOOTHINGFACTOR, self.smoothingFactorSlider),
        (ARG_CAVITIESDIAMETER, self.cavitiesDiameterSlider),
        (ARG_CAVITIESDEPTH, self.cavitiesDepthSlider),
        (ARG_SHELLTHICKNESS, self.shellThicknessSlider),
        (ARG_ITERATIONS, self.iterationsSlider),
        (ARG_SPACING, self.spacingSlider),
        (ARG_SHELLDISTANCE, self.shellDistanceSlider)
        ]:
      self.scriptedEffect.setParameter(pId, pElement.value)

    self.scriptedEffect.setParameter(ARG_CARVECAVITIES, self.cavitiesCheckBox.isChecked())
    self.scriptedEffect.setParameter(ARG_CREATESHELL, self.createShellCheckBox.isChecked())

    if self.outputSegmentationRadioButton.isChecked():
      self.scriptedEffect.setParameter(ARG_OUTPUTTYPE, OUTPUT_SEGMENTATION)
    elif self.outputModelRadioButton.isChecked():
      self.scriptedEffect.setParameter(ARG_OUTPUTTYPE, OUTPUT_MODEL)
    
    self.scriptedEffect.setParameter(ARG_OUTPUTMODELNODE, self.outputModelSelector.currentNodeID)

    self.disableOptions()
    self.cleanup()
  

  #
  # Effect specific methods (the above ones are the API methods to override)
  #

  def disableOptions(self):
    
    if self.scriptedEffect.parameter(ARG_CARVECAVITIES)=='True':
      self.cavitiesDiameterSlider.setEnabled(True)
      self.cavitiesDepthSlider.setEnabled(True)
    else:
      self.cavitiesDiameterSlider.setEnabled(False)
      self.cavitiesDepthSlider.setEnabled(False)

    if self.scriptedEffect.parameter(ARG_CREATESHELL)=='True':
      self.shellThicknessSlider.setEnabled(True)
      self.shellDistanceSlider.setEnabled(True)
    else:
      self.shellThicknessSlider.setEnabled(False)
      self.shellDistanceSlider.setEnabled(False)

    self.outputModelSelector.setEnabled(self.scriptedEffect.parameter(ARG_OUTPUTTYPE) == OUTPUT_MODEL)

  def onApply(self):

    if self.applyButton.text == 'Cancel':
      self.logic.requestCancel()
      return
    


    self.scriptedEffect.saveStateForUndo()

    self.applyButton.text = 'Cancel'
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

    #TODO: check if segment is not empty

    seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    segID = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()
    modID = self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE)
    outputMod = slicer.mrmlScene.GetNodeByID(self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE))

    kwargs = {
      ARG_CARVECAVITIES : self.scriptedEffect.parameter(ARG_CARVECAVITIES)=='True',
      ARG_CAVITIESDIAMETER : self.scriptedEffect.doubleParameter(ARG_CAVITIESDIAMETER),
      ARG_CAVITIESDEPTH : self.scriptedEffect.doubleParameter(ARG_CAVITIESDEPTH),

      ARG_CREATESHELL : self.scriptedEffect.parameter(ARG_CREATESHELL)=='True',
      ARG_SHELLTHICKNESS : self.scriptedEffect.doubleParameter(ARG_SHELLTHICKNESS),

      ARG_OUTPUTTYPE : self.scriptedEffect.parameter(ARG_OUTPUTTYPE),
      ARG_SMOOTHINGFACTOR : self.scriptedEffect.doubleParameter(ARG_SMOOTHINGFACTOR),

      ARG_ITERATIONS : self.scriptedEffect.doubleParameter(ARG_ITERATIONS),
      ARG_SPACING : self.scriptedEffect.doubleParameter(ARG_SPACING),
      ARG_SHELLDISTANCE : self.scriptedEffect.doubleParameter(ARG_SHELLDISTANCE)
    }
    
    self.logic.ApplyWrapSolidify(seg, segID, outputMod, **kwargs)

    qt.QApplication.restoreOverrideCursor()
    self.applyButton.text = 'Apply'


  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update



class WrapSolidifyLogic(object):

  def __init__(self, scriptedEffect):
    self.scriptedEffect = scriptedEffect
    self.logCallback = None
    self.cancelRequested = False
  
  def requestCancel(self):
    logging.info("User requested cancelling.")
    self.cancelRequested = True

  def ApplyWrapSolidify(self, segmentationNode, segmentID, modelNode=None, **kwargs):
    """Applies the Shrinkwrap-Raycast-Shrinkwrap Filter, a surface filter, to the selected passed segment.
    
    Arguments:
        segmentationNode (vtkMRMLSegmentationNode): Segmentation Node the filter gets applied to.
        segmentID (string): ID of the Segment the filter gets applied to. WARNING: The filter replaces the former segmentation.
        modelNode (vtkMRMLModelNode): If 'MODEL' is selected as outputType, the polydata of this node gets updated.
        **outputType (string): Possible options: 'MODEL', 'SEGMENTATION'
        **smoothingFactor (double): 0-1, Smoothing of the surface representation of the input used by this filter. This factor is also used on the output surface model and segmentation.
        **carveCavities (bool): Solidification process also carves out deeper caves.
        **cavitiesDiameter (double): 0.1-100, Entrance diameter of caves. Only used if carveCavities == True.
        **cavitiesDepth (double): 0.1-1000, Depth of caves. Only used if carveCavities == True.
        **createShell (bool): Whether or not the solidification will be done only for a outer shell.
        **shellThickness (-0.1-10): Thickness of the shell. Only used if createShell == True. WARNING: If <0, a nonmanifold mesh gets created, Segmentation will fail.
        **iterationsNr (1-10): Nr. of iterations of the shrinkwrap process.
        **spacing (>1): Spacing of remesh process.
        **shellDistance (-0.1-10): Maximal distance between input segmentation and shell, larger distant vertices get deleted. Only used if createShell == True. If <0, no vertex gets deleted.



    
    Returns:
        bool: Should return True.
    """

    self.cancelRequested = False
    self.segLogic = slicer.vtkSlicerSegmentationsModuleLogic
    self.modelsLogic = slicer.modules.models.logic()

    options = {
      ARG_OUTPUTTYPE : DEFAULT_OUTPUTTYPE,
      ARG_SMOOTHINGFACTOR : DEFAULT_SMOOTHINGFACTOR,

      ARG_CARVECAVITIES : DEFAULT_CARVECAVITIES,
      ARG_CAVITIESDIAMETER : DEFAULT_CAVITIESDIAMETER,
      ARG_CAVITIESDEPTH : DEFAULT_CAVITIESDEPTH,

      ARG_CREATESHELL : DEFAULT_CREATESHELL,
      ARG_SHELLTHICKNESS : DEFAULT_SHELLTHICKNESS,

      ARG_ITERATIONS : DEFAULT_ITERATIONS,
      ARG_SPACING : DEFAULT_SPACING,
      ARG_SHELLDISTANCE : DEFAULT_SHELLDISTANCE
      }

    options.update(kwargs)



    if self.logCallback: self.logCallback('Initializing Filtering Process...')
    
    segmentationNode.GetSegmentation().SetConversionParameter(slicer.vtkBinaryLabelmapToClosedSurfaceConversionRule().GetSmoothingFactorParameterName(), str(options[ARG_SMOOTHINGFACTOR]))
    segmentationNode.RemoveClosedSurfaceRepresentation()
    segmentationNode.CreateClosedSurfaceRepresentation()
    segmentationNode.Modified()
    if segmentationNode.GetDisplayNode():
      segmentationNode.GetDisplayNode().UpdateScene(slicer.mrmlScene)
    segmentationNode.Modified()
    segment = segmentationNode.GetSegmentation().GetSegment(segmentID)
    inputPolyData = segment.GetRepresentation(slicer.vtkSegmentationConverter().GetSegmentationClosedSurfaceRepresentationName())

    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(inputPolyData)
    cellLocator.BuildLocator()
    self.outputModel = modelNode


    #region Helper Functions

    def cleanup():
      if self.logCallback: self.logCallback('')

    
    
    def polydataToModel(polydata, smooth=True):

      if self.cancelRequested:
        cleanup()
        return False
      if self.logCallback: self.logCallback('Creating Model...')
      
      if smooth:
        polydata = smoothPolydata(polydata)


      
      if not self.outputModel:
        self.outputModel = self.modelsLogic.AddModel(polydata)
        self.outputModel.SetName(segment.GetName())
        self.outputModel.GetDisplayNode().SliceIntersectionVisibilityOn()
      

      else:
        self.outputModel.SetAndObservePolyData(polydata)
        
      self.outputModel.GetDisplayNode().SetColor(segment.GetColor())

      return True

    def polydataToSegment(polydata, smooth=True):
      if self.cancelRequested:
        cleanup()
        return False
      if self.logCallback: self.logCallback('Updating Segmentation...')
      
      if smooth:
        polydata = smoothPolydata(polydata)

      tempSegment = vtkSegmentationCorePython.vtkSegment()
      tempSegment.SetName(segment.GetName())
      tempSegment.SetColor(segment.GetColor())
      tempSegment.AddRepresentation(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(), polydata)
      segment.DeepCopy(tempSegment)
      segmentationNode.Modified()
      segmentationNode.RemoveClosedSurfaceRepresentation()
      segmentationNode.CreateClosedSurfaceRepresentation()
      
      if segmentationNode.GetDisplayNode():
        segmentationNode.GetDisplayNode().UpdateScene(slicer.mrmlScene)
      
      
      return True

    def shrinkPolydata(polydata, inputPolydata, distance=0):
      if distance == 0:
        smoothFilter = vtk.vtkSmoothPolyDataFilter()
        smoothFilter.SetInputData(0, polydata)
        smoothFilter.SetInputData(1, inputPolydata)
        smoothFilter.Update()
        return smoothFilter.GetOutput()
      
      else:
        points = polydata.GetPoints()

        for i in range(points.GetNumberOfPoints()):
          originPoint = np.array(points.GetPoint(i))
          closestPoint = np.array([0.0,0.0,0.0])
          cell = vtk.vtkGenericCell()
          cellId = vtk.mutable(0)
          subId = vtk.mutable(0)
          closestPointDist2 = vtk.mutable(0)

          cellLocator.FindClosestPoint(originPoint, closestPoint, cell, cellId, subId, closestPointDist2)

          vector = closestPoint - originPoint
          vectorLength = np.linalg.norm(vector)

          if VAL_OFFSET > 0 and vectorLength > 0.01:
            newLocation = closestPoint - ((vector/vectorLength) * VAL_OFFSET)
          else:
            newLocation = closestPoint
          
          points.SetPoint(i, newLocation)
        
        polydata.SetPoints(points)
        return polydata
    
    def remeshPolydata(polydata, spacing):
      whiteImage = vtk.vtkImageData()
      bounds = [0]*6
      polydata.GetBounds(bounds)

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
      smootherSinc.SetInputData(polydata)
      smootherSinc.SetNumberOfIterations(20)
      smootherSinc.FeatureEdgeSmoothingOff()
      smootherSinc.BoundarySmoothingOff()
      smootherSinc.NonManifoldSmoothingOn()
      smootherSinc.NormalizeCoordinatesOn()
      smootherSinc.Update()

      return smootherSinc.GetOutput()

    #endregion


    #region Initialize Shrinkwrap

    # create sphere
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

    shrinkModelPD.DeepCopy(shrinkPolydata(cleanPolyData.GetOutput(), inputPolyData, VAL_OFFSET))
    #shrinkModelPD.DeepCopy(cleanPolyData.GetOutput())

    #endregion

    
    if options[ARG_CARVECAVITIES]:
      
      #region Initialization Shrinkwrap
      
      for x in range(VAL_ITERATIONSFIRST):

        # remesh
        if self.cancelRequested:
          cleanup()
          return False
        if self.logCallback: self.logCallback('Remeshing %s/%s...' %(x+1, VAL_ITERATIONSFIRST))
        
        shrinkModelPD.DeepCopy(remeshPolydata(shrinkModelPD, [VAL_SPACINGFIRST]*3))


        # shrink
        if self.cancelRequested:
          cleanup()
          return False

        if self.logCallback: self.logCallback('Shrinking %s/%s...' %(x+1, VAL_ITERATIONSFIRST))

        shrinkModelPD.DeepCopy(shrinkPolydata(shrinkModelPD, inputPolyData, VAL_OFFSET))

      #endregion

      #region Raycast
      if self.cancelRequested:
        cleanup()
        return False
      if self.logCallback: self.logCallback('Raycasting...')

      # Find Large Faces and remember IDs of connected points
      largeCellIds = vtk.vtkIdList() # IDs of cells
      for i in range(shrinkModelPD.GetNumberOfCells()):
        cell = shrinkModelPD.GetCell(i)

        # get Length longest edge of cell
        pointsArray = list()
        for p in range(cell.GetNumberOfPoints()):
          pointsArray.append(np.array(cell.GetPoints().GetPoint(p)))

        edgeLength = list()
        for pa in range(len(pointsArray) - 1):
          length = np.linalg.norm(pointsArray[pa] - pointsArray[pa + 1])
          edgeLength.append(length)

        if max(edgeLength) > options[ARG_CAVITIESDIAMETER]:
          largeCellIds.InsertNextId(i)

      # extract large cells for cell point localization

      largeCellsPolyData = vtk.vtkPolyData()
      largeCellsPolyData.DeepCopy(shrinkModelPD)
      largeCellsPolyData.BuildLinks()

      for c in range(largeCellsPolyData.GetNumberOfCells()):
        if largeCellIds.IsId(c) == -1:
          largeCellsPolyData.DeleteCell(c)

      largeCellsPolyData.RemoveDeletedCells()

      # subdivide
      ids = vtk.vtkIdFilter()
      adapt = vtk.vtkAdaptiveSubdivisionFilter()
      adapt.SetInputData(shrinkModelPD)
      adapt.SetMaximumEdgeLength(options[ARG_CAVITIESDIAMETER]/VAL_SUBDIVISIONFRACTIONS)
      adapt.SetMaximumTriangleArea(vtk.VTK_INT_MAX)
      adapt.SetMaximumNumberOfPasses(vtk.VTK_INT_MAX)
      adapt.Update()

      clean = vtk.vtkCleanPolyData()
      clean.SetInputData(adapt.GetOutput())
      clean.Update()

      shrinkModelPD.DeepCopy(clean.GetOutput())

      if self.cancelRequested:
        cleanup()
        return False

      if largeCellIds.GetNumberOfIds() > 0 and options[ARG_CAVITIESDEPTH] > 0.0:

        # locate the points of previous large cells and write into largePointIds Set
        largeDistance = vtk.vtkImplicitPolyDataDistance()
        largeDistance.SetInput(largeCellsPolyData)

        largePointIds = set()
        for p in range(shrinkModelPD.GetNumberOfPoints()):
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

        for i in range(shrinkModelPD.GetNumberOfCells()):
          cell = shrinkModelPD.GetCell(i)
          pointIds = cell.GetPointIds()
          
          if cell.GetCellType() == 5: # cell with face
            for p in range(pointIds.GetNumberOfIds()):
              pointId = pointIds.GetId(p)

              # check if cell with point was large before subdividion, and if point got checked already
              if not pointId in largePointIds or ((pointId in vert_location_dict) and vert_location_dict[pointId][0] == False):
                
                # random False value, won't be moved
                vert_location_dict.update({pointId:(False,np.array([0.0,0.0,0.0]),0.0)})

              else:
                cell = shrinkModelPD.GetCell(i)
                pointId = cell.GetPointIds().GetId(p)
                normal = np.array(shrinkModelPD.GetPointData().GetArray('Normals').GetTuple(pointId)) * (-1)
                vector = normal * options[ARG_CAVITIESDEPTH] # max Length greater 0, checked above

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

        for i in range(numberOfPoints):
          # check result
          if vert_location_dict[i][0] == True:

            shrinkModelPD.GetPoints().SetPoint(i, vert_location_dict[i][1])
            
            # # check distance between two new locations with positive result
            # cellIds = vtk.vtkIdList()
            # shrinkModelPD.GetPointCells(i, cellIds)
            # pointChanged = False
            # for c in range(cellIds.GetNumberOfIds()):
            #   if pointChanged == True:
            #     break
            #   cell = shrinkModelPD.GetCell(cellIds.GetId(c))
            #   pointIds = cell.GetPointIds()
            #   for p in range(pointIds.GetNumberOfIds()):
            #     if pointChanged == True:
            #       break
            #     pointId = pointIds.GetId(p)
            #     if pointId != i and vert_location_dict[pointId][0] == True:
            #       point = vert_location_dict[pointId][1]
            #       distance = np.linalg.norm(-vert_location_dict[i][1] + point)
            #       if distance < options[ARG_RAYCASTMAXHITDISTANCE]:
            #         shrinkModelPD.GetPoints().SetPoint(i, vert_location_dict[i][1])
            #         pointChanged = True
      
      #endregion



    #region Main Shrinkwrap

    for x in range(int(options[ARG_ITERATIONS])+1):
      
      # remesh
      if self.cancelRequested:
        cleanup()
        return False
      if self.logCallback: self.logCallback('Remeshing %s/%s...' %(x+1, int(options[ARG_ITERATIONS]+1)))
      
      shrinkModelPD.DeepCopy(remeshPolydata(shrinkModelPD, [options[ARG_SPACING]]*3))

      if x == int(options[ARG_ITERATIONS]) and options[ARG_CREATESHELL]:
        break

      # shrink
      if self.cancelRequested:
        cleanup()
        return False
      if self.logCallback: self.logCallback('Shrinking %s/%s...' %(x+1, int(options[ARG_ITERATIONS] if options[ARG_CREATESHELL] else int(options[ARG_ITERATIONS]))))

      shrinkModelPD.DeepCopy(shrinkPolydata(shrinkModelPD, inputPolyData))

    shrinkModelPD.DeepCopy(smoothPolydata(shrinkModelPD))

    #endregion

    if options[ARG_CREATESHELL]:
      
      #region Remove Caps
      if options[ARG_SHELLDISTANCE] >= 0:
        
        if self.cancelRequested:
          cleanup()
          return False
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
          for p in range(points.GetNumberOfPoints()):
            point = points.GetPoint(p)
            distance = implicitDistance.EvaluateFunction(point)

            if abs(distance) > options[ARG_SHELLDISTANCE]:
              nonsolidPolyData.DeleteCell(c)
              break

        nonsolidPolyData.RemoveDeletedCells()
        shrinkModelPD.DeepCopy(nonsolidPolyData)

      #endregion


      #region Solidification
      if options[ARG_SHELLTHICKNESS] >= 0:
        
        if self.cancelRequested:
          cleanup()
          return False
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
          for i in range(cellIDs.GetNumberOfIds()):
            n = []
            n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3))
            n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3 + 1))
            n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3 + 2))

            normalsArray.append(np.array(n) * (-1))

          # calculate position of new vert
          dir_vec = np.zeros(3)
          
          for n in normalsArray:
            dir_vec = dir_vec + np.array(n)

          dir_vec_norm = dir_vec / np.linalg.norm(dir_vec)
          proj_length = np.dot(dir_vec_norm, np.array(normalsArray[0]))
          dir_vec_finallenght = dir_vec_norm * proj_length
          vertex_neu = np.array(shrinkModelPD.GetPoint(pointID)) + (dir_vec_finallenght * options[ARG_SHELLTHICKNESS])
          
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

        shrinkModelPD = vtk.vtkPolyData()
        shrinkModelPD.DeepCopy(triangleFilter.GetOutput())

      #endregion



    #region Create Output
    if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
      polydataToSegment(shrinkModelPD, False)
    elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
      polydataToModel(shrinkModelPD, False)
    else:
      logging.error('Unknown Output Type')
      cleanup()
      return False
      
    cleanup()
    return True

    #endregion





ARG_SMOOTHINGFACTOR = 'smoothingFactor'
DEFAULT_SMOOTHINGFACTOR = 0.2

ARG_CARVECAVITIES = 'carveCavities'
DEFAULT_CARVECAVITIES = False
ARG_CAVITIESDIAMETER = 'cavitiesDiameter'
DEFAULT_CAVITIESDIAMETER = 20.0
ARG_CAVITIESDEPTH = 'cavitiesDepth'
DEFAULT_CAVITIESDEPTH = 100.0

ARG_CREATESHELL = 'createShell'
DEFAULT_CREATESHELL = False
ARG_SHELLTHICKNESS = 'shellThickness'
DEFAULT_SHELLTHICKNESS = 1.5

ARG_OUTPUTTYPE = 'outputTypeWrapSolidify'
OUTPUT_MODEL = 'MODEL'
OUTPUT_SEGMENTATION = 'SEGMENTATION'
DEFAULT_OUTPUTTYPE = OUTPUT_SEGMENTATION

ARG_OUTPUTMODELNODE = 'outputModelNode'
DEFAULT_OUTPUTMODELNODE = ''

ARG_ITERATIONS = 'iterationsNr'
DEFAULT_ITERATIONS = 6
ARG_SPACING = 'spacing'
DEFAULT_SPACING = 1
ARG_SHELLDISTANCE = 'shellDistance'
DEFAULT_SHELLDISTANCE = 0.7


# hidden options:
VAL_ITERATIONSFIRST = 2
VAL_SPACINGFIRST = 10
VAL_OFFSET = 15
VAL_SUBDIVISIONFRACTIONS = 3


