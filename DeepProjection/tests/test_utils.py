from DeepProjection.utils import equal_float, check_img_for_nan, check_for_blank_img, check_img_for_right_shape
import numpy as np


def test_equal_float():
    """ Tests equal_float() function in utils.py """

    # Initialize variables we will be using
    a = None
    b = None

    # 1. Test to see if two "equal" floats return True.
    a = 0.0
    b = 0.0
    assert equal_float(a, b) == True

    # 2. Test to see if two different floats return False.
    a = 0.0
    b = 0.01
    assert equal_float(a, b) == False


def test_check_img_for_nan():
    """ Tests check_img_for_nan() in utils.py """

    # 1. Make a 32 x 32 numpy image with a 'nan' value in it. Function should return True.
    test_1_img = np.ones(shape=(32, 32), dtype=np.float64)
    test_1_img[6, 6] = np.nan
    assert check_img_for_nan(test_1_img) == True

    # 2. Make a 32 x 32 numpy image without any 'nan' values in it. Function should return False.
    test_2_img = np.ones(shape=(32, 32), dtype=np.float64)
    assert check_img_for_nan(test_2_img) == False


def test_check_for_blank_img():
    """ Tests check_for_blank_img() in utils.py """

    # 1. Make a 32 x 32 numpy image with 0 values in it. Function should return True.
    test_1_img = np.zeros(shape=(32, 32), dtype=np.float64)
    assert check_for_blank_img(test_1_img) == True

    # 2. Make a 32 x 32 numpy image with 1 values in it. Function should return True.
    test_2_img = np.ones(shape=(32, 32), dtype=np.float64)
    assert check_for_blank_img(test_2_img) == False


def test_check_img_for_right_shape():
    """ Tests check_img_for_right_shape() in utils.py """

    # 1. Make a 32 x 32 numpy image, pass in a shape of (32, 32). Function should return true.
    test_1_img = np.zeros(shape=(32, 32), dtype=np.float64)
    test_1_shape = (32, 32)
    assert check_img_for_right_shape(test_1_img, test_1_shape) == True

    # 2. Make a 32 x 32 numpy image, pass in a shape of (1, 32). Function should return false.
    test_2_img = np.zeros(shape=(32, 32), dtype=np.float64)
    test_2_shape = (1, 32)
    assert check_img_for_right_shape(test_2_img, test_2_shape) == False

    # 3. Make a 32 x 32 numpy image, pass in a shape of (32, 1). Function should return false.
    test_3_img = np.zeros(shape=(32, 32), dtype=np.float64)
    test_3_shape = (32, 1)
    assert check_img_for_right_shape(test_3_img, test_3_shape) == False


