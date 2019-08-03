import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

class SegmentEditorSRSFilter(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    import string
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "SegmentEditorSRSFilter"
    self.parent.categories = ["Segmentation"]
    self.parent.dependencies = ["Segmentations", "Models", "SegmentStatistics"]
    self.parent.contributors = ["Sebastian Andress (LMU Munich)"]
    self.parent.hidden = True
    self.parent.helpText = "This hidden module registers the segment editor effect"
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
      This filter was developed by Sebastian Andress (LMU University Hospital Munich, Germany, Department of General-, Trauma- and Reconstructive Surgery).
      """
    slicer.app.connect("startupCompleted()", self.registerEditorEffect)

  def registerEditorEffect(self):
    import qSlicerSegmentationsEditorEffectsPythonQt as qSlicerSegmentationsEditorEffects
    instance = qSlicerSegmentationsEditorEffects.qSlicerSegmentEditorScriptedEffect(None)
    effectFilename = os.path.join(os.path.dirname(__file__), self.__class__.__name__+'Lib/SegmentEditorEffect.py')
    instance.setPythonSource(effectFilename.replace('\\','/'))
    instance.self().register()

# class SegmentEditorSRSFilterWidget(ScriptedLoadableModuleWidget):

#   def setup(self):
#     ScriptedLoadableModuleWidget.setup(self)

class SegmentEditorSRSFilterTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_SRSFilter1()

  def test_SRSFilter1(self):
    """
    Basic automated test of the segmentation method:
    - Create segmentation by placing sphere-shaped seeds
    - Run segmentation
    - Verify results using segment statistics
    The test can be executed from SelfTests module (test name: SegmentEditorSRSFilter)
    """

    self.delayDisplay("Starting test_SRSFilter1")

    import vtkSegmentationCorePython as vtkSegmentationCore
    import vtkSlicerSegmentationsModuleLogicPython as vtkSlicerSegmentationsModuleLogic
    import SampleData
    from SegmentStatistics import SegmentStatisticsLogic

    ##################################
    self.delayDisplay("Load master volume")

    masterVolumeNode = SampleData.downloadSample('MRBrainTumor1')

    ##################################
    self.delayDisplay("Create segmentation containing a two spheres")

    segmentationNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode')
    segmentationNode.CreateDefaultDisplayNodes()
    segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(masterVolumeNode)
    
    filterModes = ["CONVEXHULL", "RAYCASTS", "DEEPHULL", "SOLIDIFIED"]
    spheres = [
      [10, 4, 4, 4],
      [10, -4,-4,-4]]
    appender = vtk.vtkAppendPolyData()
    for sphere in spheres:
      sphereSource = vtk.vtkSphereSource()
      sphereSource.SetRadius(sphere[0])
      sphereSource.SetCenter(sphere[1], sphere[2], sphere[3])
      appender.AddInputConnection(sphereSource.GetOutputPort())

    for m in filterModes:
      segmentName = str(m)
      segment = vtkSegmentationCore.vtkSegment()
      segment.SetName(segmentationNode.GetSegmentation().GenerateUniqueSegmentID(segmentName))
      appender.Update()
      segment.AddRepresentation(vtkSegmentationCore.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(), appender.GetOutput())
      segmentationNode.GetSegmentation().AddSegment(segment)

    ##################################
    self.delayDisplay("Create segment editor")

    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.show()
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
    segmentEditorNode = slicer.vtkMRMLSegmentEditorNode()
    slicer.mrmlScene.AddNode(segmentEditorNode)
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setMasterVolumeNode(masterVolumeNode)

    ##################################
    self.delayDisplay("Run SRS-Filter Effect")
    segmentEditorWidget.setActiveEffectByName("SRS-Filter")
    effect = segmentEditorWidget.activeEffect()

    for t in ["MODEL", "SEGMENTATION"]:
      effect.setParameter("outputType", t)
      for m in filterModes:
        self.delayDisplay("Creating Output Type %s, Filter Mode %s" %(t,m))
        effect.setParameter("filterMode", m)
        segmentEditorWidget.setCurrentSegmentID(segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(m))
        effect.self().onApply()

    ##################################
    self.delayDisplay("Creating Segments from Models")
    for m in filterModes:
      model = slicer.util.getNode(m)
      segmentName = "MODEL_%s" % m
      segment = vtkSegmentationCore.vtkSegment()
      segment.SetName(segmentationNode.GetSegmentation().GenerateUniqueSegmentID(segmentName))
      segment.SetColor(model.GetDisplayNode().GetColor())
      segment.AddRepresentation(vtkSegmentationCore.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(), model.GetPolyData())
      segmentationNode.GetSegmentation().AddSegment(segment)

    ##################################
    self.delayDisplay("Compute statistics")
    segStatLogic = SegmentStatisticsLogic()
    segStatLogic.getParameterNode().SetParameter("Segmentation", segmentationNode.GetID())
    segStatLogic.getParameterNode().SetParameter("ScalarVolume", masterVolumeNode.GetID())
    segStatLogic.computeStatistics()
    statistics = segStatLogic.getStatistics()

    ##################################
    self.delayDisplay("Check a few numerical results")
    
    self.assertEqual( round(statistics["CONVEXHULL",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 47896)
    self.assertEqual( round(statistics["MODEL_CONVEXHULL",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 47896)

    self.assertEqual( round(statistics["RAYCASTS",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 41628)
    self.assertEqual( round(statistics["MODEL_RAYCASTS",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 41628)

    self.assertEqual( round(statistics["DEEPHULL",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 6450)
    self.assertEqual( round(statistics["MODEL_DEEPHULL",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 6450)

    self.assertEqual( round(statistics["SOLIDIFIED",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 2477)
    self.assertEqual( round(statistics["MODEL_SOLIDIFIED",'ScalarVolumeSegmentStatisticsPlugin.volume_mm3']), 2477)
    
    self.delayDisplay('test_SRSFilter1 passed')
