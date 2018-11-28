import numpy as np

def extract_label(proatlas_data,label,threshold=0):
    '''
    extract special label map from atlas
    Parameters:
    proatlas : 4d array,1-3 dimension repersents xyz coorodinate,4th dimension repersents different label
    label    :  extract map with corresponding label
    threshold:  extract map which have been deleted values under threshold value.
    '''
    coordinate_shape = proatlas_data.shape[:3]
    extract_atlas = np.zeros(coordinate_shape)
    extract_atlas[proatlas_data[:,:,:,label] > threshold ] = 1
    return extract_atlas

def extract_threshold(promap_data,threshold=0):
    '''
    extract special map data below a certain threshold.
    Parameters:
    promap   :  3d array,1-3 repersents xyz coorodinates
    threshold:  extract map which will be deleted values below threshold value.
    '''
    coordinate_shape = promap_data.shape
    extract_map = np.zeros(coordinate_shape)
    extract_map[promap_data[:, :, :] > threshold] = 1
    return extract_map

def extract_sp_val(promap_data,sp_val):
    '''extract special value map:
    Parameters:
    promap    :  3d array,1-3 repersents xyz coorodinates
    sp_val :  special value will be extracted and binaryzed.'''
    coordinate_shape = promap_data.shape
    extract_map = np.zeros(coordinate_shape)
    for i in sp_val:
        extract_map[promap_data[:, :, :] == i] = 1
    return extract_map

