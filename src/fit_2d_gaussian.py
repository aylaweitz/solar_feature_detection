import sys
sys.path.insert(1, '/Users/weitz/Documents/LMSAL-BAERI-work-repo')

from imports import *
from load_data_functions import *
from dot_detection_functions import *


####### 2D GAUSSIAN FUNCTION ########
def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = xy
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))

    return g.ravel()


####### FITTING 2D GAUSSIAN ########
def fit_2d_gaussian(one_dot, frame, dataset, R, plotting=False, bounds=None):
    """
    one_dot:

    frame: 

    dataset

    R

    plotting: bool, default False (optional)
        plot Gaussian fit.

    bounds: 2d list, default (optional)
    """

    xcenter, ycenter = one_dot[frame]

    data = dataset[frame].data

    data_subset = np.array(data[ycenter-R:ycenter+R, xcenter-R:xcenter+R]) # square subset, size 2R

    # making meshgrid? idk
    x = np.linspace(-R, R, 2*R)
    y = np.linspace(-R, R, 2*R)
    xx,yy = np.meshgrid(x, y)

    # Flatten the meshgrid and data_subset to fit with curve_fit expectations
    x_flat = xx.ravel()
    y_flat = yy.ravel()
    data_flat = data_subset.ravel()

    # initial guess
    initial_guess = [5, 0, 0, 1, 1, 10, 1]

    if bounds == None:
        bounds = ([-10, -2, -2, -2, -2, -100, -10],
                [10, 2, 2, 2, 2, 100, 10])

    try:
        popt, pcov = scipy.optimize.curve_fit(twoD_Gaussian, (x_flat, y_flat), data_flat, p0=initial_guess, bounds=bounds)
        print(popt)
    except:
        if plotting == True:
            plt.imshow(data_subset, cmap=plt.cm.jet, origin='lower',
                extent=(x.min(), x.max(), y.min(), y.max()))
            plt.show()
        print('could not fit :(')
        return np.nan, np.nan, False

    data_fitted = twoD_Gaussian((xx, yy), *popt)

    if plotting == True:
        plt.imshow(data_subset, cmap=plt.cm.jet, origin='lower',
                extent=(x.min(), x.max(), y.min(), y.max()))
        plt.contour(x, y, data_fitted.reshape(data_subset.shape), 5, colors='w')
        plt.show()

    return popt, pcov, True


# REPORT KM SIZE
# GET SIZES
def get_size(one_dot, peak_i, dataset, R, plotting=False, bounds=None):

    popt, pcov, is_success = fit_2d_gaussian(one_dot, peak_i, dataset, R, plotting, bounds)

    if is_success:

        # get x and y sizes
        sigma_x = popt[3]
        sigma_y = popt[4]

        # rotation angle
        theta = popt[5]

        # Calculate potential standard deviations along the rotated axes
        sigma_1 = np.sqrt(sigma_x**2 * np.cos(theta)**2 + sigma_y**2 * np.sin(theta)**2)
        sigma_2 = np.sqrt(sigma_x**2 * np.sin(theta)**2 + sigma_y**2 * np.cos(theta)**2)

        # Identify major and minor axes
        sigma_major = max(sigma_1, sigma_2)
        sigma_minor = min(sigma_1, sigma_2)

        # Convert to FWHM
        FWHM_major = 2 * np.sqrt(2 * np.log(2)) * sigma_major
        FWHM_minor = 2 * np.sqrt(2 * np.log(2)) * sigma_minor

        # errors 
        err = np.diag(pcov)
        sigma_x_err = err[3]
        sigma_y_err = err[4]

        sigma_err = np.sqrt(sigma_x_err**2 + sigma_y_err**2)

        # FWHM
        # fwhm_x = 2 * np.sqrt(2 * np.log(2)) * sigma_x
        # fwhm_y = 2 * np.sqrt(2 * np.log(2)) * sigma_y
    
        # convert to km
        km_size_major = convert_pix_to_dist(dataset[0], FWHM_major)
        km_size_minor = convert_pix_to_dist(dataset[0], FWHM_minor)

        return km_size_major, km_size_minor, sigma_err

    else:
        return np.nan, np.nan, np.nan
    