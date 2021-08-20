import sys
import numpy as np
import skopi as sk

import skimage.measure as sm
from numba import jit


############################################ From calculate_diffraction_image_resolution.ipynb ############################################

def calculate_maximum_diffraction_resolution(beam_file, pdb_file, detector_dimensions, increase_factor=1):
    """
    Calculates the maximum diffraction resolution for a set of PDB images.
    
    The calculations in this function are done for images generated with skopi's SimpleSquareDetector only.
    
    Parameters
    ----------
    beam_file: str
        Directory path to beam file used to simulate diffraction images whose resolution we are trying
        to calculate.
    pdb_file: str
        Directory path to PDB file used to simulate diffraction images whose resolution we are trying to calculate.
    detector_dimensions: tuple(int, float, float)
        Tuple representing the dimensions of SimpleSquareDetector used to simulate diffraction images whose
        resolution we are trying to calculate.
        Tuple format is (num pixels for row and col, detector size, detector distance from protein).
    increase_factor: int (optional)
        Multiplies the number of photons/pulse by given factor. Default value does not affect photons/pulse.
        
    Return
    ------
    max_resolution: float
        Maximum resolution possible for diffraction images simulated with beam, PDB, and detector dimensions.
    """
    
    ## Maximum resolution to return.
    max_resolution = None

    # Define the detector dimensions
    n_pixels, det_size, det_dist = detector_dimensions
    
    # Setup only the beam, particle, and detector, as they contain the information
    # we need to calculate the scatter angles
    beam = sk.Beam(beam_file)
    if increase_factor != 1:
        beam.set_photons_per_pulse(increase_factor * beam.get_photons_per_pulse())
    particle = sk.Particle()
    particle.read_pdb(pdb_file, ff = 'WK')
    det = sk.SimpleSquareDetector(int(n_pixels), float(det_size), float(det_dist), beam=beam)
    
    # Calculate the maximum resolution for a diffraction image
    pixel_reciprocal_vectors = np.squeeze(det.pixel_distance_reciprocal)
    max_resolution = 1.0 / pixel_reciprocal_vectors[0][0]
    
    # Return our result
    return max_resolution

############################################### From check_for_nan_values.ipynb #########################################################


# From check_for_nan_values.ipynb
def equal_float(a, b):
    """ Helper function to see if two floats are equal. """
    return abs(a - b) <= sys.float_info.epsilon


# From check_for_nan_values.ipynb
def check_img_for_nan(img):
    """
    Checks for 'nan' values in an image.
    
    Parameters
    ----------
    img: numpy.array
        Image to check for 'nan'
        
    Returns
    -------
    True, if img does have a 'nan'.
    False, if img does not have a 'nan'.
    """

    if np.isnan(img).any():
        return True
    else:
        return False


# From check_for_nan_values.ipynb
def check_for_blank_img(img):
    """
    Checks to see if the given image is blank.
    
    Parameters
    ----------
    img: numpy.array
        Image to check to see if it is blank.
    
    Return
    ------
    True, if image is blank.
    False, if image is not blank.
    """

    # Take absolute of img before taking the sum in order to
    # prevent corner case of negative and positive values resulting in
    # a sum of zero.
    if equal_float(np.sum(np.absolute(img)), 0):
        return True
    else:
        return False


# From check_for_nan_values.ipynb
def check_img_for_right_shape(img, shape_to_check):
    """
    Checks to see if given image has the right shape.
    
    Parameters
    ----------
    img: numpy.array
        Image to check shape.
    shape_to_check: tuple
        Shape to check against given image's shape.
    
    Returns
    -------
    True, if image has the specified shape.
    False, if image does not have the specified shape.
    """
    
    if img.shape == shape_to_check:
        return True
    else:
        return False
            
############################################ From modified_thumbnail_generation.ipynb #####################################################
# Upsampling
@jit(nopython=True)
def upsample(warr, dim, binr, binc):
    """ Helper function to img_to_thumbnail() """
    upCalib = np.zeros(dim)
    #for k in range(dim[0]):
    for i, ix in enumerate(range(0,dim[0],binr)):
        if ix+binr > dim[0]:
            er = dim[0]+1
        else:
            er = ix+binr
        for j, jy in enumerate(range(0,dim[1],binc)):
            if jy+binc > dim[1]:
                ec = dim[1]+1
            else:
                ec = jy+binc
            upCalib[ix:er,jy:ec] = warr[i,j]
    return upCalib


# perform binning here
def downsample(assem, bin_row=2, bin_col=2, mask=None):
    """ Helper function to img_to_thumbnail() """
    if mask is None:
        combinedMask = np.ones_like(assem)
    else:
        combinedMask = mask
    downCalib = sm.block_reduce(assem, block_size=(bin_row, bin_col), func=np.sum)
    downWeight = sm.block_reduce(combinedMask, block_size=(bin_row, bin_col), func=np.sum)
    warr = np.zeros_like(downCalib, dtype='float32')
    ind = np.where(downWeight > 0)
    warr[ind] = downCalib[ind] / downWeight[ind]
    return warr


def img_to_thumbnail(img, beam_file, pdb_file, detector_dimensions, n_part_per_shot, thumbnail_shape):
    """
    Generates an image into a thumbnail image.
    
    Works only for images generated using skopi's SimpleSquareDetector.
    
    Parameters
    ----------
    img: numpy.array
        Image, with shape of (row, col), to create a thumbnail out of.
    beam_file: str
        Path to beam file used to image img.
    pdb_file: str
        Path to PDB file used to make img.
    detector_dimensions: tuple(int, float, float)
        Tuple representing the dimensions of SimpleSquareDetector used to simulate diffraction images whose
        resolution we are trying to calculate.
        Tuple format is (num pixels for row and col, detector size, detector distance from protein).
    n_part_per_shot: int
        Number of particles that were simulated in the image.
    thumbnail_shape: tuple(int, int)
        Tuple representing the dimension of the thumbnail image.
        Tuple format is (thumbnail_rows, thumbnail_cols).
    
    Return
    ------
    thumbnail_img: numpy.array
        Thumbnail version of image passed in.
    """
    
    # Thumbnail image to return.
    thumbnail_img = None
    
    # Set increase factor, as seen in original thumbnail generation notebook.
    increase_factor = 1000
    
    # Get detector dimension variables.
    n_pixels, det_size, det_dist = detector_dimensions
    
    # Set the number of rows and columns for binning.
    # Formula for rows is: bin_rows = img.shape[0] / thumbnail_shape[0], or bin_rows = (rows in img) / (desired thumbnail rows).
    # Formula for cols is: bin_cols = img.shape[1] / thumbnail_hape[0], or bin_cols = (cols in img) / (desired thumbnail rows).
    bin_rows = img.shape[0] / thumbnail_shape[0]
    bin_cols = img.shape[1] / thumbnail_shape[1]
    
    # Setup beam file.
    beam = sk.Beam(beam_file)
    beam.set_photons_per_pulse(increase_factor * beam.get_photons_per_pulse())

    # Setup particle file.
    particle = sk.Particle()
    particle.read_pdb(pdb_file, ff='WK')
    
    # Setup detector.
    det = sk.SimpleSquareDetector(int(n_pixels), float(det_size), float(det_dist), beam=beam)

    # Setup experiment.
    exp = sk.SPIExperiment(det, beam, particle, n_part_per_shot)

    # Get a pattern from experiment. This is used to create a mask to make thumbnail out of.
    pattern = exp.generate_image_stack(return_photons=True, return_intensities=False)
    
    # Get mask to be used for making the thumbnail.
    mask = det.assemble_image_stack(np.ones_like(pattern))
    
    # Generate the thumbnail.
    thumbnail_img = downsample(det.assemble_image_stack(img), bin_row=bin_rows, bin_col=bin_cols, mask=mask)
    
    # Return thumbnail image.
    return thumbnail_img
