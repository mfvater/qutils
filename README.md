# QuTILs
Based on the publication 'QuTILs: open-source image-based infiltrating immune cell detection for research application' at: https://doi.org/10.1038/s41523-026-01021-z

Standalone workflow is provided to produce raw values used in downstream calculations for total TILs %, estimated stromal invasive marginal TILs %, and hotspot clusters. Classifier and associated files can be incorporated in other workflows/pipelines as desired. 

## Pre-Processing

1) Create a new QuPath project, placing the 'classifiers' folder with object_classifiers and density_maps sub-folders, as well as the 'scripts' folder into the directory. This allows the software to appropriately find the necessary components to run the workflow. Then load target images into the project.
2) It is a good idea to first calculate the estimated stain vectors for the image(s) in the project. Foregoing this step can influence the detection, feature calculation, classification, and downstream annotation accuracy/consistency.
3) Either manually segment the tissue ROI or generate a whole-image annotation region. In the latter case, this will add significant processing time with the inclusion of brightfield background, and may also include undesirable regions (smudging, dye discoloration, trapped air, sharpie marks etc.) that will likely unfortunately show up as detected/classified objects. This can be refined by usage of a trained pixel classifier, not provided.
4) Workflow script can then be run and results will append to a .txt file in the main project directory. Processing speed will depend on local vs. external drive storage of image files, and users on very limited hardware may encounter memory-based errors depending on the size of the tissue ROI and image file.

## Citation

If you use this tool in your research, please cite the parent publication.

Vater, M., Salgado, R., Blige, E., Lefebvre, H., Strahan, A., Gatti-Mays, M., Rugo, H.S., Ballman, K., Watson, M., Wen, Y. and Leon-Ferre, R., 2026. QuTILs: open-source image-based infiltrating immune cell detection for research application. npj Breast Cancer. https://doi.org/10.1038/s41523-026-01021-z
