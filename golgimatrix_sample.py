
#Copyright, Dmitry Yampolsky 2022
#golgimatrix project

#code sample

#Kai M. Bracey, Kung-Hsien Ho, Dmitry Yampolsky, Guogiang Gu, Irina Kaverina, William R. Holmes,
#Microtubules Regulate Localization and Availability of Insulin Granules in Pancreatic Beta Cells,
#Biophysical Journal, Volume 118, Issue 1, 2020, Pages 193-206

#dependencies:

# CircStat2012a
# SPHERE_VORONOI
# resize
# GridSphere
# UniformSampling
# imgaussian

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter
from skimage.transform import resize
from skimage.filters import threshold_otsu
from tkinter import filedialog
import os

# Add necessary paths
import sys
sys.path.append('./CircStat2012a')
sys.path.append('./SPHERE_VORONOI')
sys.path.append('./resize')
sys.path.append('./GridSphere')
sys.path.append('./UniformSampling')

def golgimatrixwin12():
    noised_switch = False
    generate_m_switch = True
    ncases = 1
    matrix_size = 2**2
    nmatrixes = 5

    draw_switch = False

    npixels = 1000
    bin_ceil = 314

    # Load data
    filename = filedialog.askopenfilename(filetypes=[("GIF files", "*.gif")])
    if not filename:
        return

    cdata = plt.imread(filename)
    cdata = np.squeeze(cdata)

    if cdata.ndim != 3:
        raise ValueError('Image file must be a stack')

    if len(np.unique(cdata)) != 2:
        cdata = cdata.astype(float)
        cdata2 = resize(cdata, (100, 100, cdata.shape[2]))

        thresh_seek = 2
        thresh_mark = -1

        pix_diff = -3
        pix_diff2 = -3

        while abs(pix_diff) > 2:
            while np.sign(pix_diff2) == np.sign(pix_diff):
                pix_diff = pix_diff2
                thresh_seek += thresh_mark
                cdata2 = irish_threshold(cdata, 100, thresh_seek)
                pix_diff2 = np.sum(cdata2) - npixels

            thresh_mark = -thresh_mark * 0.5
            pix_diff = pix_diff2

        cdata = cdata2
    else:
        cdata = cdata / np.max(cdata)  # make binary

        if np.prod(cdata.shape) / 2 < np.sum(cdata):  # if background (majority) is black
            cdata = ~cdata  # invert

    # Image data is now binary

    if cdata.shape[0] != cdata.shape[1]:  # make image size square
        cdata = make_square(cdata)

    cdata_original = cdata

    # Trans to vector
    three_d_matrix_coords = make_spherical_vectors(cdata)

    results_array = []

    uniform_model = particle_sample_sphere(40)
    iod, d_hist = iod_from_norm(three_d_matrix_coords, uniform_model)

    max_freq = 50
    d_hist_ceiled = np.minimum(d_hist, max_freq)

    print('IOD = ')
    print(np.var(d_hist_ceiled) / np.mean(d_hist_ceiled))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(three_d_matrix_coords[:, 0], three_d_matrix_coords[:, 1], three_d_matrix_coords[:, 2], c='b', marker='.')
    plt.show()

def make_spherical_vectors(binary_in):
    x, y, z = np.meshgrid(range(binary_in.shape[0]), range(binary_in.shape[1]), range(binary_in.shape[2]))

    # 0:1 -> -1:1
    x = 2 * (x / binary_in.shape[0]) - 1
    y = 2 * (y / binary_in.shape[1]) - 1
    z = 2 * (z / binary_in.shape[2]) - 1

    # Make 3d coord array
    vec_tmp = np.column_stack((x[binary_in], y[binary_in], z[binary_in]))
    # Normalize to unit sphere
    vec_out = vec_tmp / np.sqrt(np.sum(vec_tmp**2, axis=1))[:, np.newaxis]
    return vec_out

def make_square(data_in):
    max_dim = max(data_in.shape[:2])
    diff_tmp = abs(data_in.shape[0] - data_in.shape[1])
    diff_tmp1 = int(np.ceil(diff_tmp / 2))
    diff_tmp2 = int(np.floor(diff_tmp / 2))

    data_out = np.zeros((max_dim, max_dim, data_in.shape[2]))
    if data_in.shape[0] < data_in.shape[1]:
        data_out[diff_tmp1:-diff_tmp2, :, :] = data_in
    else:
        data_out[:, diff_tmp1:-diff_tmp2, :] = data_in
    return data_out

def geodesic_grid(scale_in, rand_factor_in):
    # This function needs to be implemented based on the GridSphere function
    pass

def make_gauss_model(n_in, res_in, sigma_in, r_in, thresh_in):
    uniform_model_gauss = particle_sample_sphere(n_in)
    uniform_model_gauss = (uniform_model_gauss + 1) / 2
    uniform_model_gauss = np.floor(uniform_model_gauss * res_in).astype(int)
    zero_mesh = np.zeros((res_in, res_in, res_in))

    for i in range(n_in):
        zero_mesh[tuple(uniform_model_gauss[i])] = 1

    uniform_model_gauss2 = gaussian_filter(zero_mesh, sigma=sigma_in, truncate=r_in)
    uniform_model_gauss2_bin = (uniform_model_gauss2 > thresh_in)
    vec_out = make_spherical_vectors(uniform_model_gauss2_bin)
    geo = geodesic_grid(len(vec_out), 0)
    return vec_out

def random_distribution_model(npts):
    ex = 2 * np.random.rand(npts, 3) - 1  # -1:1 cube of random points
    mask = np.sum(ex**2, axis=1) < 1
    ex2 = ex[mask]  # cutout the unit sphere
    pts_out = ex2 / np.sqrt(np.sum(ex2**2, axis=1))[:, np.newaxis]  # move to sphere surface
    return pts_out

# Note: Some functions like irish_threshold, particle_sample_sphere, and iod_from_norm
# are not defined here and would need to be implemented separately.

if __name__ == "__main__":
    golgimatrixwin12()