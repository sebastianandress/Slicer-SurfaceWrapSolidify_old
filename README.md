**Copyright (c) 2019, Sebastian Andreß**\
All rights reserved. Please find the license [here](https://github.com/sebastianandress/Slicer-SegmentEditorSRSFilter/blob/master/LICENSE.md).

For further collaborations, patient studies or any help, do not hesitate to contact [Sebastian Andreß, MD](sebastian.andress@med.uni-muenchen.de).

# Shrinkwrap-Raycast-Shrinkwrap (SRS) Surface Filter

## Introduction
This filter was designed for creating fractured bone models for fast 3D printing. Especially in orthopedic trauma surgery, the editing time, as well as the printing time should be as short as possible. Using this filter helps to fulfil both features. Also, by removing inner spongious structures, it is possible to archive a fracture reduction on the printed model.

In our use-case, we used this filter after applying a simple threshold operation and separating the bone with simple brushing and island techniques. Please watch the [workflow example](#Workflow-Example) videos. The filter was tested on more than 30 acetabular fracture models, it reduced the printing time about 70%.

![](/Resources/Screenshots/screenshot1.png)

## Principle Description and Parameters
The SRS-Filter uses the following pipeline to archive the filter result. Also the related parameters are explained.

1. A surface representation of the selected segment is created (segmented model).
    * __Smoothing Factor__: Smoothing of the surface representation. This factor is also used on the output surface model and segmentation.

2. Around this segmented model a sphere is created (sphere/surface model).

3. Iteratively shrinkwrapping and remeshing the sphere model to the segmented model, keeping a distance between both models.

    ![](/Resources/Media/description_1_firstSR.png)

    * __Offset First Shrinkwrap__: Distance between both models.\
    _space between red and turquoise lines/segmentation_
    * __Iterations First Shrinkwrap__: The more iterations, the more precise this first result will be. _WARNING_: Will increase processing time, depending on spacing.
    * __Spacing First Remesh__: The resolution of the sphere model. The smaller the spacing and higher the resolution, the more precise this first output result will be. _WARNING_: Will drastically increase processing time. Especially for the first iteration, high resolution is not needed in most cases.\
    _Space between white vertex points_

4. Since holes deeper then broad will not be reached by shrinkwrapping, they are covered with long-edged faces. These faces are detected and separated. The resulting vertices are raycasted along their normals to the hitting point on the segmented model, if they fulfil specific criterias described in the following parameters.

    ![](/Resources/Media/description_2_raycastexplanation.png)
    ![](/Resources/Media/description_2_raycast.png)

    * __Raycast Search Edge Length__: If one edge of a face is longer than this, it is used in this step.\
    _Image of step 3: length of red line between the two white points above the joint-cup representation_
    * __Raycast Output Edge Length__: Length of these edges after subdivision. The shorter this length, the more the raycasting will possibly hit spongious parts inside fracture gaps.\
    _Upper left: length of red lines between white points_
    * __Raycast Max. Hit Distance__: Maximal distance between two hits. This setting is to prevent accidental hits way of the other hits.\
    _Upper right: grey circles around hits_
    * __Raycast Max./Min. Length__: Maximal/Minimal length of the ray.\
    _Upper right: length of blue dotted line_


5.  After an initial remeshing operation, another shrinkwrap-rememesh iteration finally approximates the former sphere model to the segmented model.

    ![](/Resources/Media/description_3_secondSR.png)

    * __Iterations Second Shrinkwrap__: The more iterations, the more precise this final result will be. _WARNING_: Will increase processing time, depending on spacing.
    * __Spacing Second Remesh__: The resolution of the surface model. The smaller the spacing and higher the resolution, the more precise the resulting surface model will be. _WARNING_: Will drastically increase processing time. A value of 1 mm^3 normally is enough depending on the original (CT/MRI) volume spacing.\
    _Space between white vertex points_

6. All vertices and faces of the surface model, that finally are not close to the segmented model, are deleted, resulting in a nonmanifold model.

    ![](/Resources/Media/description_4_nonmanifold.png)

    * __Max. Models Distance__: If a vertex of the surface model is farer away from the segmented model, it gets deleted with it faces.
    
7. A solidifying operation creates a solid, manifold model.

    ![](/Resources/Media/description_5_solidified.png)

    * __Solidification Thickness__: Thickness of the solidified surface model. The solidification is performed. It is performed to the inside of the model.

8. Optional:  A segmentation is created stencilling the solidified model.

### Kwargs for code implementation
| Parameter | Key | Default | Possible Values | Unit | Type |
| - | - | - | - | - | - |
| Offset First Shrinkwrap | `offsetFirstShrinkwrap` | `15.0` | `>0` | mm | `float` |
| Spacing First Remesh | `spacingFirstRemesh` | `10.0` | `>0` | mm^3 | `float` |
| Spacing Second Remesh | `spacingSecondRemesh` | `1.0` | `>0` | mm^3 | `float` |
| Iterations First Shrinkwrap | `iterationsFirstShrinkwrap` | `3` | `>0` |  |  `int` |
| Iterations Second Shrinkwrap | `iterationsSecondShrinkwrap` | `5` | `>=0` |  | `int` |
| Raycast Search Edge Length | `raycastSearchEdgeLength` | `20.0` | `>0` | mm | `float` |
| Raycast Output Edge Length | `raycastOutputEdgeLength` | `2.0` | `>0` | mm | `float` |
| Raycast Max. Hit Distance | `raycastMaxHitDistance` | `2.0` | `>0` | mm | `float` |
| Raycast Max. Length | `raycastMaxLength` | `100.0` | `>0` | mm | `float` |
| Raycast Min. Length | `raycastMinLength` | `0.0` | `>=0` | mm | `float` |
| Max. Models Distance | `maxModelDistance` | `0.7` | `>=0` | mm | `float` |
| Solidification Thickness | `solidificationThickness` | `1.5` | `>0` | mm | `float` |
| Smoothing Factor | `smoothingFactor` | `0.2` | `0-1` | % | `float` |
| Output Type | `outputType` | `MODEL` | `MODEL`, `SEGMENTATION` |  | `string` |
| Filter Mode | `filterMode` | `SOLIDIFIED` | `CONVEXHULL`, `RAYCASTS`, `DEEPHULL`, `NONMANIFOLD`, `SOLIDIFIED` |  | `string` |


## Results
To see the creation of the results, please see the [Workflow Example](#Workflow-Example) section.

After a common threshold procedure and slightly manual edits, the resulting segmentation shows parts of the cartilage and spongious bone.
 
![](/Resources/Media/result_seg_threshold.gif)

Applying the SRS-Filter with the default settings results in the following segmentation.

![](/Resources/Media/result_seg_srs.gif)

The following interim results are possible:

### Convex Hull
![](/Resources/Media/result_mod_convexhull.gif)

### Raycasts
![](/Resources/Media/result_mod_raycast.gif)

### Deep Hull
![](/Resources/Media/result_mod_deephull.gif)

### Nonsolid Model
![](/Resources/Media/result_mod_nonsolid.gif)

### Solidified Surface
![](/Resources/Media/result_mod_solid.gif)

## Workflow Example

As described in the [Introduction](#Introduction), this filter was designed to create fast printable bones as easy and least time consuming as possible. The parameters are especially fitted for hemipelvic bones. The threshold and manual edit takes about 3 minutes, the filter itself another 2 minutes resulting in an printable bone. The filter itself reduces the printing time about 70 percent.

### Threshold
![](/Resources/Media/workflow_threshold.gif)

A thresholding operation between 300 and the maximal Houndsfiled units was performed, using the __Threshold__ Effect. By using the sphere brush, first the femoral head, and subsequently connecting parts in the sacroiliac joint were erased. Using the __Islands__ Effect, the exempted hemipelvis was added to an own segment.

### Filter process
![](/Resources/Media/workflow_processing.gif)

In this example, the processing time was 1:46 min on a Apple MacBook Pro 2017 (3,1 GHz Intel Core i7, Memory 16 GB).


## How to install
Please follow the description on the official [3D Slicer page](https://www.slicer.org/wiki/Documentation/Nightly/Developers/FAQ/Extensions). The extension should also be soon available in the builtin Extension Manager.