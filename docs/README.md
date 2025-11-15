# Documentation

This folder contains educational reference material for understanding the theory and implementation of stellar cluster photometry and density profile analysis.

## Theory and Methods

### [BACKGROUND_ESTIMATION_THEORY.md](BACKGROUND_ESTIMATION_THEORY.md)
Comprehensive guide to CCD background estimation for photometry:
- Physical sources of background (zodiacal light, airglow, scattered light)
- Why spatial variations matter
- 2D mesh-based sigma-clipping algorithm (detailed math)
- Impact on photometric uncertainties
- Validation examples
- Comparison with alternative methods

**Use this for**: Understanding how background subtraction works, implementing your own background estimation, troubleshooting photometry issues.

### [PIPELINE_WORKFLOW.md](PIPELINE_WORKFLOW.md)
End-to-end workflow for stellar cluster analysis:
- Architecture: modular .py files + visualization notebooks
- Complete data flow: FITS → photometry → membership → density profile
- Integration of all analysis steps
- File structure and organization
- How each notebook connects to the next

**Use this for**: Understanding the overall analysis structure, seeing how the pieces fit together, planning your own cluster analysis project.

## How to Use These Resources

These documents are designed to be:
1. **Educational** - Learn the theory and best practices
2. **Reference** - Look up equations, parameters, and methods
3. **Implementation guides** - See how theory translates to code

They complement the Jupyter notebooks by providing deeper theoretical context and explaining design decisions.

## Related Notebooks

- `data_reduction_simple.ipynb` - Implements the workflow described in PIPELINE_WORKFLOW.md
- `membership.ipynb` - Uses methods described in membership sections of both documents
- `completeness.ipynb` - Background estimation impacts completeness modeling

## Questions or Improvements?

If you find errors or have suggestions for improving these documents, please open an issue or submit a pull request!
