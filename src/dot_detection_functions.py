import sys
sys.path.insert(1, '/Users/weitz/Documents/LMSAL-BAERI-work-repo')

from imports import *
from load_data_functions import *


# MAKING UNSHARP MASK FUNCTION
def unsharp_mask(data, blur_r):
    """
    Perform unsharp masking on image -- gaussian filter.
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.filters.gaussian_filter.html

    Parameters
    ----------
    data: np.array
        data to enhance
    blur_r: int
        radius to blur by

    Returns
    -------
    np.array
        unsharp masked image
    """
    blurred = scipy.ndimage.filters.gaussian_filter(data, sigma=blur_r, mode='constant')
    result = data - blurred

    # https://scikit-image.org/docs/stable/auto_examples/filters/plot_unsharp_mask.html
    # result = skimage.filters.unsharp_mask(data, radius=blur_r, amount=amount)
    # im2 = im1.filter(ImageFilter.UnsharpMask(radius = blur_r, percent = 200))

    return result


def make_unsharp_masked_maps(obs, blur_r):
    """
    Make array of unsharp masked data.

    Parameters
    ----------
    obs: np.array of SunPy maps

    blur_r: int
        radius to blur by

    Returns
    -------
    np.array of SunPy maps
        unsharp masked maps (header info is kept the same)
    """
    masked_maps = [sunpy.map.Map(unsharp_mask(one_map.data, blur_r), one_map.wcs) for one_map in obs]
    masked_maps = sunpy.map.Map(masked_maps, sequence=True)
    return masked_maps


def calc_fwhm(dist, gauss):
    # calculate the full-width half max of a gaussian
    try:
        spline = UnivariateSpline(dist, gauss - np.max(gauss)/2, s=0)
        r1, r2 = spline.roots() # find the roots -- x positions where gauss is at half of its maximum value
        return r2 - r1 # subtract the two locations to get full width half max
    except:
        return np.nan



class dot_obs():
    """one observation of a dot"""
    def __init__(self, dot, time, earth_time, frame, image, thresh, instrument):
        self.y = dot[0]
        self.x = dot[1]
        self.r = dot[2]
        self.time = time
        self.earth_time = earth_time
        self.frame = frame
        self.image = image
        self.thresh = thresh # the threshold this dot observation was found at (e.g. 3 stds)
        self.instrument = instrument

    def __str__(self):
        return f'({self.x}, {self.y}), rad: {self.r}, frame: {self.frame}, time: {self.time}'
        
        
    def check_if_same_dot(self, other_dot_obs, dist=2):
        """Check if other dot observation is of the same dot

        Parameters
        ----------
        other_dot_obs: dot_obs object
            other dot observation to check
        dist: integer, optional (default 2)
            max distance (in pixels) away to be considered the same dot

        Returns
        -------
        True if other blob is less than dist away from current blob
        False otherwise
        """
        dist_to_blob = np.sqrt( (self.y - other_dot_obs.y)**2 + (self.x - other_dot_obs.x)**2 )
        
        # check time too?
        
        if dist_to_blob < dist: # if distance to other blob is less than 2 pix -- arbitrary val?? maybe make it less than or equal to???
            return True, dist_to_blob
        else:
            return False, dist_to_blob

    
    def search_other_wavelength(self, other_data, other_time_data, quietx, quiety, search_r=2, thresh=2, blur_r=2, instrument=''):
        """Search around dot to see if visible in other wavelength
        dist: distance from center of dot in initial wavelength detection"""

        # create search region centered around dot observation
        searchx = [self.x-search_r, self.x+search_r] 
        searchy = [self.y-search_r, self.y+search_r]

        other_wavelength_dots = search_for_dots(self.frame, other_data, other_time_data, searchx, searchy, quietx, quiety, thresh, blur_r, instrument)

        return other_wavelength_dots

        
    def convert_pix_to_dist(self, pix_val):
        """
        Convert pixel value to distance
        """
        pixscale_x = self.image.scale[0] # deg/pix in x
        pixscale_y = self.image.scale[1] # deg/pix in y
        d_sun = self.image.dsun # Observer distance from the center of the Sun
        return ((pix_val*u.pix * pixscale_x).to(u.rad) * d_sun) / (u.rad)
    
    
    def fit_one_gaussian(self, r_add=3, plotting=False):
        # get intensity over horizontal and vertical cuts

        if np.isnan(self.x): # if no dot detection, cannot fit gaussian so set variables to np.nan
            self.popt_h = np.nan
            self.popt_v = np.nan
            self.pcov_h = np.nan
            self.pcov_v = np.nan
            self.fwhm_h = np.nan
            self.fwhm_v = np.nan
            self.enhancement = np.nan
            return

        else: # if these is a dot detection, fit gaussian
        
            data = np.array(self.image.data)

            r = self.r + r_add # detected radius of the dot plus number of pixels specified by r_add

            x, y = int(self.x), int(self.y)
            x1, x2 = int(self.x-r), int(self.x+r) + 1
            y1, y2 = int(self.y-r), int(self.y+r) + 1

            horizontal_cut = data[y, x1:x2]
            peak_i = np.argmax(horizontal_cut) # center horizontal peak (assume vertical peak is in similar location)
            middle_i = len(horizontal_cut)//2 # get middle index
            diff_i = peak_i - middle_i # difference in index values between peak and middle
            x_start = x1 + diff_i; y_start = y1 + diff_i # add in the difference in index values to center peak -- can only center with integer precision
            x_end = x2 + diff_i + 1; y_end = y2 + diff_i + 1 # +1 to include end point

            # make new cuts with peak intensity centered
            horizontal_cut = data[y, x_start:x_end]
            vertical_cut = data[y_start:y_end, x]

            # normalize by max value
            horizontal_cut = horizontal_cut / np.max(horizontal_cut)
            vertical_cut = vertical_cut / np.max(vertical_cut)

            # converted pixels to km!
            dist = [self.convert_pix_to_dist(i).to(u.km).value for i in range(len(horizontal_cut))]

            # FIT GAUSSIANS TO DATA -- asymmetic model, sum of two gaussians
            # def gaussian2(x, amp1, mean1, stddev1, amp2, mean2, stddev2, cst):
            #     return amp1 * np.exp(-((x - mean1) / 4 / stddev1)**2) + amp2 * np.exp(-((x - mean2) / 4 / stddev2)**2) + cst
            
            def gaussian1(x, amp, mean, stddev):
                return amp * np.exp(-((x - mean) / 4 / stddev)**2)

            def exp_gaussian(x, amp, mean, stddev, amp2, f, g, cst): # non-constant background gaussian -- gaussian + exponential
                return amp * np.exp(-((x - mean) / 4 / stddev)**2) + amp2 * f**(-x/g) + cst

            # initial conditions
            p0 = [0.8, # amp
                1300, # mean
                100, # stddev
                1, # amp2
                1.02, # f
                200, # g
                0 # cst
                ]

            # set bounds for params --- amp, mean, stddev, amp2, f, g, cst
            bounds = ([0.1, 50, 50, -2, -2, 1, 0], [1, 1500, 100, 3, 3, 300, 0.5])

            try: # try to fit gaussian
                horizontal_popt, horizontal_pcov = scipy.optimize.curve_fit(exp_gaussian,  # horizontal
                                                                xdata=dist, 
                                                                ydata=horizontal_cut, 
                                                                p0 = p0, bounds = bounds)
                horizontal_gaussian = exp_gaussian(np.array(dist), *horizontal_popt)

                vertical_popt, vertical_pcov = scipy.optimize.curve_fit(exp_gaussian, # vertical
                                                                            xdata=dist, 
                                                                            ydata=vertical_cut, 
                                                                            p0 = p0, bounds = bounds)
                vertical_gaussian = exp_gaussian(np.array(dist), *vertical_popt)

            except: # if can't fit gaussian, set vals to np.nan
                print('cant fit gaussian')
                horizontal_popt, horizontal_pcov = np.nan, np.nan
                vertical_popt, vertical_pcov = np.nan, np.nan
                
                horizontal_gaussian = [np.nan for val in dist]
                vertical_gaussian = [np.nan for val in dist]
                
                one_horizontal_gaussian = [np.nan for val in dist]
                one_vertical_gaussian = [np.nan for val in dist]
                
                
            # CALCULATE FWHM
            # sigma_h = horizontal_popt[2] # standard deviation of distribution
            horizontal_fwhm = calc_fwhm(dist, horizontal_gaussian)#-np.min(horizontal_gaussian))
            
            # sigma_v= vertical_popt[2] # standard deviation of distribution
            vertical_fwhm = calc_fwhm(dist, vertical_gaussian)#-np.min(vertical_gaussian))
        

            # ESTIMATE INTENSITY ENHANCEMENT -- done on gaussian fit or the data itself??
            # get two ends of the intensity profile
            horizontal_peak_intensity = np.max(horizontal_cut)
            horizontal_background = np.mean([horizontal_cut[0], horizontal_cut[-1]])
            horizontal_enhancement = horizontal_peak_intensity / horizontal_background # unsure is correct formula

            vertical_peak_intensity = np.max(vertical_cut)
            vertical_background = np.mean([vertical_cut[0], vertical_cut[-1]])
            vertical_enhancement = vertical_peak_intensity / vertical_background

            avg_enhancement = np.mean([horizontal_enhancement, vertical_enhancement])


            # make gaussian plot
            # plotting
            if plotting == True:
                fig = plt.figure()
                plt.title(f'Average Intensity Enhancement {avg_enhancement:.3}')
                plt.xlabel('Distance [km]')
                plt.ylabel('Normalized Intensity')
                plt.scatter(dist, horizontal_cut, 
                            c='r', label=f'Horizontal Cut\nFWHM = {horizontal_fwhm:.5} km')
                plt.plot(dist, horizontal_gaussian, c='r', zorder=0, alpha=0.5)

                plt.scatter(dist, vertical_cut,
                            c='k', label=f'Vertical Cut\nFWHM = {vertical_fwhm:.5} km')
                plt.plot(dist, vertical_gaussian, c='k', zorder=0, alpha=0.5)
                
                plt.xlim(dist[0], dist[-1])
                plt.legend()
                plt.show()


            # set parameters for this dot_obs
            self.popt_h = horizontal_popt
            self.popt_v = vertical_popt
            self.pcov_h = horizontal_pcov
            self.pcov_v = vertical_pcov
            self.fwhm_h = horizontal_fwhm
            self.fwhm_v = vertical_fwhm
            self.enhancement = avg_enhancement
            self.fitting_r = r_add # the addition tail used for gaussian fitting
            # self.gaussian_plot = fig
         
        
    def fit_gaussian_multiple_attempts(self, rs = [2, 3, 4, 5, 6, 7], plotting=False):

        lowest_uncerts = []
        lowest_cov = []
        params_to_keep = []
        lowest_uncert_fwhm = []
        lowest_r = np.nan

        match = False
        
        for r in rs: # try each additional pixels to consider for gaussian fitting, see which one does the best fit
            # print(r)
            self.fit_one_gaussian(r_add=r, plotting=plotting)
            
            if self.popt_h is not np.nan: # if there was a successful fit

                match = True # there's a match, so set match=True

                param_uncerts = [ [np.sqrt(self.pcov_h[i,i]) for i in range(len(self.pcov_h))], [np.sqrt(self.pcov_v[i,i]) for i in range(len(self.pcov_v)) ] ]
                covs = [self.pcov_h, self.pcov_v]
                params = [self.popt_h, self.popt_v]
                fwhms = [self.fwhm_h, self.fwhm_v]

                if len(lowest_uncerts) != 0: # if there has already been a successful fit

                    # check if the uncertainty of these parameters are smaller or larger than the previously determined best fit, (this fit - previous lowest error fit)
                    uncerts_change_h = [param_uncerts[0][elem_i] - lowest_uncerts[0][elem_i] for elem_i in range(len(param_uncerts[0]))] # horizontal
                    uncerts_change_v = [param_uncerts[1][elem_i] - lowest_uncerts[1][elem_i] for elem_i in range(len(param_uncerts[1]))] # vertical

                    # if mean of uncerts_change_h/v is negative, that means (previous lowest error fit) > (this fit) so WE ARE intered in keeping this new fit
                    if np.nanmean([uncerts_change_h, uncerts_change_v]) < 0:
                        # print(f'better fit: {np.nanmean([uncerts_change_h, uncerts_change_v])}')
                        lowest_uncerts = param_uncerts
                        lowest_covs = covs
                        params_to_keep = params
                        lowest_uncert_fwhm = fwhms
                        lowest_r = r
                    # else, do not change the lowest uncertainty parameters

                else: # if there has not alreayd been a successful fit, set this fit as the lowest
                    # print('first fit, so initializing')
                    lowest_uncerts = param_uncerts
                    lowest_covs = covs
                    params_to_keep = params
                    lowest_uncert_fwhm = fwhms
                    lowest_r = r

        if match == True:
            # once finished looping through all r's to try, set the lowest uncertainty parameters
            self.popt_h = params_to_keep[0]
            self.popt_v = params_to_keep[1]
            self.pcov_h = lowest_uncerts[0]
            self.pcov_v = lowest_uncerts[1]
            self.fwhm_h = lowest_uncert_fwhm[0]
            self.fwhm_v = lowest_uncert_fwhm[1]
            self.fitting_r = lowest_r
        
                
#                 # calculate relative uncertainties on fit parameters
#                 uncert_h = np.sqrt(np.diag(self.pcov_h))
#                 rel_uncert_h = np.array([uncert_h[i]/self.popt_h[i] for i in range(len(uncert_h))])
                
#                 uncert_v = np.sqrt(np.diag(self.pcov_v))
#                 rel_uncert_v = np.array([uncert_v[i]/self.popt_v[i] for i in range(len(uncert_v))])
                
#                 rel_uncert = np.concatenate((rel_uncert_h, rel_uncert_v))
                
#                 # if all relative uncertainties are under 200%
#                 if np.where(abs(rel_uncert) > 2)[0].size == 0:
#                     print(f'Good enough fit! Average relative uncertainty of {np.nanmean(np.abs(rel_uncert))*100:.2f}%')
#                     break
#                 else:
#                     print(f'Too high uncertainties of fit params: {rel_uncert}')
#                     continue
                
            # print('Not able to fit gaussian')


        
    
    
class dot():
    """observations of one dot groups together"""
    def __init__(self, dot_obs_list):
        self.dot_obs_list = dot_obs_list
        
    def get_xs(self):
        return np.array([one_dot.x for one_dot in self.dot_obs_list])
    
    def get_ys(self):
        return np.array([one_dot.y for one_dot in self.dot_obs_list])
    
    def get_times(self):
        return np.array([one_dot.time for one_dot in self.dot_obs_list])

    def get_earth_times(self):
        return np.array([one_dot.earth_time for one_dot in self.dot_obs_list])

    def get_vel(self):
        xs = self.get_xs()
        x_diffs = np.diff(xs) * u.km # CONVERT TO KM

        ys = self.get_ys()
        y_diffs = np.diff(ys) * u.km

        times = self.get_times()
        time_diffs = np.diff(times) / np.timedelta64(1, 's') * u.s # s

        vel_x = x_diffs / time_diffs
        vel_y = y_diffs / time_diffs
        return [vel_x.to(u.km/u.s), vel_y.to(u.km/u.s)] # km / s


    def get_avg_vel(self):

        xs = self.get_xs()
        ys = self.get_ys()
        times = self.get_times()

        # get nan indexes if there are
        if np.isnan(xs).any():
            nan_i = np.where(np.isnan(xs))[0]

            xs = np.delete(xs, nan_i) # get rid of nan entries
            ys = np.delete(ys, nan_i)
            times = np.delete(times, nan_i)

            if len(xs) < 2:
                return [np.nan, np.nan] # km /s


        time_diff = (times[-1] - times[0]) / np.timedelta64(1, 's') * u.s # seconds

        if time_diff > 0: # ensure there is at least 2 observations at different times

            dist_pix_x = xs[-1] - xs[0] # pixels
            dist_pix_y = ys[-1] - ys[0] # pixels
            dist_x = self.dot_obs_list[0].convert_pix_to_dist(dist_pix_x).to(u.km) # km
            dist_y = self.dot_obs_list[0].convert_pix_to_dist(dist_pix_y).to(u.km) # km

            vel_x = dist_x / time_diff
            vel_y = dist_y / time_diff
            return [vel_x.value, vel_y.value] # km / s

        else:
            return [np.nan, np.nan]


    def get_frames(self):
        return [one_dot.frame for one_dot in self.dot_obs_list]

    def get_thresh(self):
        return [one_dot.thresh for one_dot in self.dot_obs_list]

    def get_peak_i(self):
        # get index of peak intensity observation -- need to run fit_gaussians first
        peak_i = np.argmax(self.get_enhancement())
        return peak_i

        
    def add_obs(self, dot_obs, position='end'):
        # add dot observation to the list of dots

        if position == 'end': # default, add dot observation to the exisiting list of dot observations
            index = len(self.dot_obs_list)

        elif position == 'start': # add observation to the beginning
            index = 0

        (self.dot_obs_list).insert(index, dot_obs) # insert new dot observation at the index specified -- if nothing given, put observation at the end
        

    def get_lifetime(self):
        time_start = self.dot_obs_list[0].time
        time_end = self.dot_obs_list[-1].time
        
        return (time_end - time_start).to(u.s)


    def search_around_detected_dot(self, when, dist, obs, time_obs, quietx, quiety, thresh, blur_r):
    
        if when == 'start': # checking frames BEFORE
            index = 0
            frame_change = -1
            
        elif when == 'end': # checking frames AFTER
            index = -1
            frame_change = +1
        
        unmatched_frames = 0
        while unmatched_frames < 2: # no more than 1 frame where the dot does not appear

            # search frame before dot appears around first location of the dot
            first_x, first_y = self.get_xs()[index], self.get_ys()[index]
            first_frame = self.get_frames()[index]
            next_frame = first_frame + frame_change # either one before first frame (-1) or one after last frame (+1)

            if next_frame < 0 or next_frame > len(obs)-1: # if reach the beginning or end of the observation, stop searching
                break

            searchx = [int(first_x-dist), int(first_x+dist)] # create search region centered around dot observation
            searchy = [int(first_y-dist), int(first_y+dist)]

            before_detection = search_for_dots(frame = next_frame, 
                                                obs = obs,
                                                time_obs = time_obs,
                                                searchx = searchx,
                                                searchy = searchy,
                                                quietx = quietx, # consider same quiet region as before
                                                quiety = quiety,
                                                thresh = thresh,
                                                blur_r = blur_r) # same blurring as previously


            if len(before_detection) > 0: # if at least one dot was detected

                new_obs_list = self.dot_obs_list.copy() # MAKE COPY OTHERWISE CHANGING new_obs_list WILL CHANGE THE ORIGINAL dot_obs_list

                for detect in before_detection: # in case there is more than one dot detection in this region in this frame
                    new_obs_list.append(detect) # append all found ones to new_obs_list

                # try grouping dot observations including the new dot observations
                updated_dots = grouping_dots(new_obs_list, dist=dist)
                the_updated_dot = 0 # initialize the_updated_dot as 0

                # if more than one dot is grouped, choose the dot with frames shared with original dot
                for a_dot in updated_dots:
                    if np.all(self.dot_obs_list) in a_dot.dot_obs_list: #np.all(a_dot.dot_obs_list) in self.dot_obs_list:
                        the_updated_dot = a_dot # choose the dot that contains the original dot_obs_list

                    else: # if does not contain original dot_obs_list
                        print(a_dot.dot_obs_list, self.dot_obs_list)
                        continue # keep the_updated_dot as np.nan

                # if the_updated_dot is not 0 (that's what it's initialized to be), check with observation to add and then add it
                if the_updated_dot != 0:

                    # if this new obs is in the group
                    matched = False
                    for a_obs in the_updated_dot.dot_obs_list: # loop through obs in the_updated_dot

                        if a_obs in before_detection: # if one of the observations is new
                            self.dot_obs_list = the_updated_dot.dot_obs_list # set the current dot's dot_obs_list to the_updated_dot's dot_obs_list
                            matched = True # set match to True

                            # only one new obs can be added (since that is how the grouping_dots function works -- one obs per frame), so can break this for loop
                            break

                    if matched: # if match was found
                        unmatched_frame = 0 # reset unmatched counter

                    else: # if a new observation is not in the group, don't add obs and unmatched_frames += 1
                        unmatched_frames += 1
                        # print('there is no new obs in group')


                else: # if when grouping the dot observations, there is no grouping that included the original grouped observations
                    unmatched_frames += 1
                    print('no group w original obs')

            else: # if not even one dot observation was detected in the next frame
                unmatched_frames += 1
                # print('no detections in before/after frame')


    def search_other_wavelength(self, other_data, other_time_data, quietx, quiety, dist=2, search_r=2, thresh=2, blur_r=2, instrument=''):
        """Search around dot to see if visible in other wavelength
        search_r: distance from center of dot in initial wavelength detection"""

        # make list of dot observations
        other_wavelength_obs_list = []

        # loop through all observation of this dot
        for one_dot_obs in self.dot_obs_list:

            frame = one_dot_obs.frame

            # search for dot in other wavelength for one observation
            other_wavelength_obs = one_dot_obs.search_other_wavelength(other_data, other_time_data, quietx, quiety, search_r, thresh, blur_r, instrument)

            if len(other_wavelength_obs) == 1: # if one detection
                other_wavelength_obs_list.append(other_wavelength_obs[0])

            elif len(other_wavelength_obs) == 0: # no detection
                continue

            else: # more than one detection, choose the one that's closest to the other wavelength obs
                i_to_keep = np.argmin([ np.sqrt((one_dot_obs.x - other_obs.x)**2 + (one_dot_obs.y - other_obs.y)**2) for other_obs in other_wavelength_obs ])
                other_wavelength_obs_list.append(other_wavelength_obs[i_to_keep])

        
        if len(other_wavelength_obs_list) > 2: # if at least 2 observations -- must have at least 2 frames that overlap 

            # GROUP DOT OBSERVATIONS
            other_wavelength_dots = grouping_dots(other_wavelength_obs_list, dist=dist)

            empty_dot = [np.nan, np.nan, np.nan] # make empty dot detection
            empty_obs = dot_obs(empty_dot, other_time_data[frame].date, get_earth_time(other_time_data[frame]), frame, other_data[frame], thresh, instrument)
            other_wavelength_dot = dot([empty_obs]) # initialize other_wavelength_dot as an empty observation

            # check if more than one dot grouped
            if len(other_wavelength_dots) > 1:
                
                # choose longest lived dot -- ???
                grouped_dot_lifetimes = [a_dot.get_lifetime().value for a_dot in other_wavelength_dots]

                # set other_wavelength_dot to this dot
                other_wavelength_dot = other_wavelength_dots[np.argmax(grouped_dot_lifetimes)]

            elif len(other_wavelength_dots) == 1: # if only one dot grouped
                other_wavelength_dot = other_wavelength_dots[0]

            # SEARCH BEFORE AND AFTER
            other_wavelength_dot.search_around_detected_dot('start', dist, other_data, other_time_data, quietx, quiety, thresh, blur_r) # search in frames BEFORE
            other_wavelength_dot.search_around_detected_dot('end', dist, other_data, other_time_data, quietx, quiety, thresh, blur_r) # search in frames AFTER

            # other_wavelength_obs_list = other_wavelength_dot.dot_obs_list # get dot obs list from the test dot object


        # if len(other_wavelength_dot) > 1: # if there is more than one matched dot
        #     dists = [ np.sqrt( (np.nanmean(self.get_xs()) - np.nanmean(one_dot.get_xs()))**2 + (np.nanmean(self.get_ys()) - np.nanmean(one_dot.get_ys()))**2 ) for one_dot in other_wavelength_dot]# choose the closest one to the other obs
        #     other_wavelength_dot = other_wavelength_dot[np.argmin(dists)]

        # elif len(other_wavelength_dot) == 1:
        #     other_wavelength_dot = other_wavelength_dot[0]

        else: # no dot
            empty_dot = [np.nan, np.nan, np.nan] # make empty dot detection

            #(self, dot, time, frame, image, thresh, instrument)
            empty_obs = dot_obs(empty_dot, other_time_data[frame].date, get_earth_time(other_time_data[frame]), frame, other_data[frame], thresh, instrument)
            other_wavelength_dot = dot([empty_obs])

        return other_wavelength_dot
    

    def get_enhancement(self):
        # intensity at peak / intensity of background
        return [one_dot.enhancement for one_dot in self.dot_obs_list]
    
    def get_fwhm(self):
        # [horizontal, vertical]
        return np.array([[one_dot.fwhm_h, one_dot.fwhm_v] for one_dot in self.dot_obs_list])

    def get_fwhm_when_peak(self):
        # [horizontal, vertical]
        fwhm = self.get_fwhm()
        peak_i = self.get_peak_i()
        peak_fwhm_h = fwhm[:, 0][peak_i]
        peak_fwhm_v = fwhm[:, 1][peak_i]
        return [peak_fwhm_h, peak_fwhm_v] # fwhm reported when dot is at peak brightness in (Alpert et al. 2016)
    
    def fit_gaussians(self, plotting=False):
        for obs in self.dot_obs_list:
            obs.fit_gaussian_multiple_attempts(plotting=plotting)

    def summary_plot(self, figsize=(12, 8)):
        fig = plt.figure(figsize=figsize)
        grid = plt.GridSpec(3, 7, wspace=0.5, hspace=0.1)

        self.fit_gaussians()

        times = self.get_times()

        fwhm_h = np.array(self.get_fwhm()[:, 0])
        fwhm_peak_h = self.get_fwhm_when_peak()[0]

        fwhm_v = np.array(self.get_fwhm()[:, 1])
        fwhm_peak_v = self.get_fwhm_when_peak()[1]

        x_pos = self.get_xs()

        intensity = self.get_enhancement()
        max_i = np.argmax(intensity)

        one_to_plot = self.dot_obs_list[max_i] # plot when dot is brightest
        peak_time = times[max_i]
        image = one_to_plot.image
        region_r = 5

        # plotting
        ax1 = fig.add_subplot(grid[:, :3], projection=one_to_plot.image.wcs)

        ax2 = fig.add_subplot(grid[0, 3], projection=one_to_plot.image.wcs)
        ax3 = fig.add_subplot(grid[1, 3], projection=one_to_plot.image.wcs)
        ax4 = fig.add_subplot(grid[2, 3], projection=one_to_plot.image.wcs)

        ax5 = fig.add_subplot(grid[0, 4:])
        ax6 = fig.add_subplot(grid[1, 4:], sharex=ax5)
        ax7 = fig.add_subplot(grid[2, 4:], sharex=ax5)

        # where is the dot located
        x1, x2 = one_to_plot.x-region_r, one_to_plot.x+region_r
        y1, y2 = one_to_plot.y-region_r, one_to_plot.y+region_r
        ax1.set_ylim(0, image.data.shape[0])
        ax1.imshow(image.data, cmap=image.cmap, norm=image.plot_settings['norm'])
        rect = patches.Rectangle((x1, y1), 
                        (x2-x1), (y2-y1),
                        color='red', linewidth=1, fill=False)
        ax1.add_patch(rect)


        # plot start of dot obs
        start_obs = self.dot_obs_list[0]
        ax2.imshow(unsharp_mask(start_obs.image.data, 2), cmap=start_obs.image.cmap, norm=plt.Normalize(-50, 100))
        ax2.set_xlim(x1, x2)
        ax2.set_ylim(y1, y2)
        ax2.set_title(f'Start\n{start_obs.time.datetime.time().replace(microsecond=0)}')
        ax2.axis('off')

        # plot peak of dot obs
        peak_obs = one_to_plot
        ax3.imshow(unsharp_mask(peak_obs.image.data, 2), cmap=peak_obs.image.cmap, norm=plt.Normalize(-50, 100))
        ax3.set_xlim(x1, x2)
        ax3.set_ylim(y1, y2)
        ax3.set_title(f'Peak\n{peak_obs.time.datetime.time().replace(microsecond=0)}')
        ax3.axis('off')

        # plot end of dot obs
        end_obs = self.dot_obs_list[-1]
        ax4.imshow(unsharp_mask(end_obs.image.data, 2), cmap=end_obs.image.cmap, norm=plt.Normalize(-50, 100))
        ax4.set_xlim(x1, x2)
        ax4.set_ylim(y1, y2)
        ax4.set_title(f'End\n{end_obs.time.datetime.time().replace(microsecond=0)}')
        ax4.axis('off')

        # plot light curve
        ax5.plot(times, intensity/np.nanmax(intensity), c='k', alpha=0.7)
        ax5.scatter(times, intensity/np.nanmax(intensity), c='k')
        ax5.axvline(peak_time, c='k', linestyle=':', zorder=0)
        ax5.set_ylabel('Intensity Compared to Background')

        # plot size over time
        ax6.scatter(times, fwhm_h, label=f'FWHM x = {fwhm_peak_h:.2f}', c='r')
        ax6.scatter(times, fwhm_v, label=f'FWHM y = {fwhm_peak_v:.2f}', c='k')
        ax6.axvline(peak_time, c='k', linestyle=':', zorder=0)
        ax6.legend()

        # plot x position over time
        ax7.scatter(times, x_pos, label=f'FWHM x = {fwhm_peak_h:.2f}', c='k')
        ax7.axvline(peak_time, c='k', linestyle=':', zorder=0)

        # times = self.get_times()

        # fwhm_h = np.array(self.get_fwhm()[:, 0])
        # fwhm_peak_h = self.get_fwhm_when_peak()[0]

        # fwhm_v = np.array(self.get_fwhm()[:, 1])
        # fwhm_peak_v = self.get_fwhm_when_peak()[1]

        # intensity = self.get_enhancement()
        # max_i = np.argmax(intensity)

        # one_to_plot = self.dot_obs_list[max_i] # plot when dot is brightest
        # image = one_to_plot.image
        # region_r = 5

        # fig, ax = plt.subplots(2, 2, figsize=(12, 6))

        # # plot dot size
        # ax[0,0].plot(times, fwhm_h, label=f'horizontal FWHM = {fwhm_peak_h:.2f}', c='r')
        # ax[0,0].axhline(fwhm_peak_h, c='r', linestyle=":")
        # ax[0,0].plot(times, fwhm_v, label=f'vertical FWHM = {fwhm_peak_v:.2f}', c='k')
        # ax[0,0].axhline(fwhm_peak_v, c='k', linestyle=":")
        # ax[0,0].legend()
        # ax[0,0].set_ylabel('Size [km]')

        # # what to plot here?
        # ax[1,0].plot(times, fwhm_h/fwhm_v, c='k', label='hor FWHM / vert FWHM')
        # ax[1,0].axhline(1, c='lightgrey', zorder=0)
        # ax[1,0].legend()
        # ax[1,0].set_ylabel('Size Ratio')

        # # plot intensity enhancement -- light curve
        # ax[0, 1].plot(times, intensity, c='k')
        # ax[0, 1].set_ylabel('Intensity Compared to Background')

        # # plot image zoomed in on dot
        # x1, x2 = one_to_plot.x-region_r, one_to_plot.x+region_r
        # y1, y2 = one_to_plot.y-region_r, one_to_plot.y+region_r
        # # ax[1,1].set_xlim(x1, x2)
        # ax[1, 1].set_ylim(0, image.data.shape[0])
        # ax[1, 1].imshow(image.data, cmap=image.cmap)
        # rect = patches.Rectangle((x1, y1), 
        #                  (x2-x1), (y2-y1),
        #                  color='red', linewidth=1, fill=False)
        # ax[1, 1].add_patch(rect)

        # ax[1,0].set_xlabel('Time')
        # # ax[1,1].set_xlabel('Time')
        # plt.show()


    def hmi_hist(self, hmi_aligned, hmi_matched, plot=False, region_r=6):
        frames = self.get_frames()
        blob_x = np.mean(self.get_xs())
        blob_y = np.mean(self.get_ys())

        # blob region
        x1, x2 = int(blob_x-region_r), int(blob_x+region_r)
        y1, y2 = int(blob_y-region_r), int(blob_y+region_r)

        hmi_time_prev = datetime.datetime.now() # set arbitrary datetime val
        all_n, all_bins = [], [] # keep list of all bin values and the counts of each bin
        for frame in frames:
            
            hmi_time = hmi_matched[frame].date.datetime.time()
            
            if hmi_time != hmi_time_prev: # if new hmi data (some repeats because EUI has faster cadence than HMI)
            
                one_image = hmi_aligned[frame]
                one_data = np.array(one_image.data)
                one_data = one_data[y1:y2, x1:x2] # cut data to region of interest

                bin_edges = np.arange(-1000, 1000, 100)
                n, bins = np.histogram(one_data, bins=bin_edges) # 10 bins

                if plot == True: # if the plot argument is set to true, plot the histogram
                    plt.hist(one_data.flatten(), bins=bins, alpha=0.3, label=hmi_time)
                    plt.xlim(-1000, 1000)
                    plt.axvline(0, linestyle='--', c='r')
                    plt.ylabel('Number of Pixels Per Bin')
                    plt.xlabel('B_z [Gauss]')
                    plt.xlim(-1000, 1000)
                    plt.legend()

                all_n.extend(n)
                all_bins.extend( [(bins[i]-bins[i+1]) + bins[i] for i in range(len(bins)-1)] )

                hmi_time_prev = hmi_time

        # checking for mixed polarity flux
        neg_flux, pos_flux = 0, 0
        for i, one_bin in enumerate(all_bins):
            if one_bin < 0: # negative flux
                neg_flux += all_n[i]
            elif one_bin > 0:
                pos_flux += all_n[i]

        if neg_flux > 0 and pos_flux > 0: # if negative and positive flux, mixed polarity
            self.magnetic_flux = 'mixed'
        elif neg_flux > 0: # if negative flux
            self.magnetic_flux = 'negative'
        elif pos_flux > 0: # if positive flux
            self.magnetic_flux = 'positive'



    # look at one dot for the entire time its detected and one frame before
    def plot(self, other_data=[], other_data_times=[], figsize=(13, 4), fontsize=5):
        
        region_r = 6 # HOW LARGE TO MAKE SURROUNDING AREA TO COMPARE TO
        
        frames = self.get_frames()

        frames_plotting = frames.copy()
        # frames_plotting.insert(0, frames[0]-plot_before_after) # plot one before
        # frames_plotting.append(frames[-1]+plot_before_after) # plot one after

        # make frames list an array
        frames = np.array(frames)

        blob_x = np.nanmean(self.get_xs())
        blob_y = np.nanmean(self.get_ys())
        
        # blob region
        x1, x2, y1, y2 = int(blob_x-region_r), int(blob_x+region_r), int(blob_y-region_r), int(blob_y+region_r)
        
        plot_size = len(frames_plotting)

        fig, axs = plt.subplots(figsize=figsize)
        axs.get_xaxis().set_visible(False); axs.get_yaxis().set_visible(False)
        axs.spines["top"].set_visible(False); axs.spines["right"].set_visible(False)
        axs.spines["bottom"].set_visible(False); axs.spines["left"].set_visible(False)

        datasets_to_include = len(other_data) + 1

        for i in range(len(frames_plotting)):
            
            frame = frames_plotting[i]

            if frame in frames:

                blob_present = True
                blob_i = np.where(frames == frame)[0][0]

                blob_x = self.get_xs()[blob_i]
                blob_y = self.get_ys()[blob_i]
                thresh = self.get_thresh()[blob_i]

            else:
                blob_present = False
            
            eui_plotting_i = i + 1 
            
            # eui
            one_image = self.dot_obs_list[blob_i].image
            one_data = np.array(one_image.data)

            ax = fig.add_subplot(datasets_to_include, plot_size, eui_plotting_i, projection=one_image.wcs)
            ax.imshow(unsharp_mask(one_data, 2), cmap=one_image.cmap, norm = plt.Normalize(0, 100))
            ax.set_xlim(x1, x2)
            ax.set_ylim(y1, y2)
            ax.set_title(self.dot_obs_list[blob_i].time.datetime.time().replace(microsecond=0), fontsize=fontsize) # get rid of microseconds in title
            ax.axis('off')

            if blob_present:
            
                if thresh == 3: # if dot found over 3 signma threshold, make marker 'lime'
                    color = 'r'
                elif thresh == 2: # if dot found over 2 sigma threshold, make marker 'red
                    color = 'lime'
                else: # if other threshold, set to black -- make these colors a normalized color map in the future?
                    color = 'k'

                ax.plot(blob_x, blob_y, marker='+', c=color)


            for num, data in enumerate(other_data):

                other_plotting_i = (num+1) * len(frames_plotting) + (i+1)

                time = other_data_times[num][frame].date.datetime.time().replace(microsecond=0)

                one_image = data[frame]
                one_data = np.array(one_image.data)

                ax = fig.add_subplot(datasets_to_include, plot_size, other_plotting_i, projection=one_image.wcs)
                ax.imshow(one_data, cmap=one_image.cmap, norm=one_image.plot_settings['norm'])
                ax.set_xlim(x1, x2)
                ax.set_ylim(y1, y2)
                ax.set_title(time, fontsize=fontsize)
                ax.axis('off')

                if blob_present:
                
                    ax.plot(blob_x, blob_y, marker='+', c=color)


        plt.tight_layout(pad=6)
        plt.show()

        return fig


# class matched_dot():
#     """observations of one dot groups together in multiple wavelengths"""
#     def __init__(self, dots):
#         self.dot_obs_list = dot_obs_list
            
        
           
    
#####


def search_for_dots(frame, obs, time_obs, searchx, searchy, quietx, quiety, thresh=2, blur_r=2, instrument=''):
    image = obs[frame]
    date = time_obs[frame].date
    earth_time = get_earth_time(time_obs[frame])
    
    # unsharp masking, blur_r x blur_r
    enhanced = unsharp_mask(image.data, blur_r)
    
    # search for dots are the base of the plume -- constrain blob search
    search_x1, search_x2 = searchx
    search_y1, search_y2 = searchy

    # to make code run faster, cut "enhanced" array so just around search area  -- include buffer in case dots on edge
    buffer = 0
    enhanced_cut = np.array(enhanced)[search_y1-buffer:search_y2+buffer, search_x1-buffer:search_x2+buffer]

    # quiet region
    quiet_x1, quiet_x2 = quietx
    quiet_y1, quiet_y2 = quiety
    subset = enhanced[quiet_y1:quiet_y2, quiet_x1:quiet_x2]
    quiet_mean = np.nanmean(subset)
    quiet_std = np.nanstd(subset)

    detection_thresh = thresh * (quiet_std) # detection threshold
    # detection_thresh = quiet_mean + thresh * quiet_std

    enhanced_cut = enhanced_cut - quiet_mean # detects bright dots on dark background -- subtract by mean of quiet region??

    blobs = blob_dog(enhanced_cut, threshold=detection_thresh) # detects bright dots on dark background -- subtract by mean of quiet region??
    blobs[:, 2] = blobs[:, 2] * np.sqrt(2)

    blobs_keep = []
    for blob in blobs:

        y = blob[0] + (search_y1-buffer) # put in terms of the entire image frame
        x = blob[1] + (search_x1-buffer)
        r = blob[2]

        true_blob = [y, x, r]

        # if y > search_y1 and y < search_y2 and x > search_x1 and x < search_x2: # if dots within search region, include them in blob_keep array
        blobs_keep.append( dot_obs(true_blob, date, earth_time, frame, image, thresh, instrument) )
        
        # else:
        #     continue
            
    return blobs_keep


def search_range_for_dots(obs, times_obs, indexes, searchx, searchy, quietx, quiety, thresh=2, blur_r=2, instrument=''):
    # TRACKING DOTS
    blobs = []
    for i in tqdm(indexes):
        one_image_dots = search_for_dots(i, obs, times_obs, searchx, searchy, quietx, quiety, thresh, blur_r, instrument)
        
        for dot in one_image_dots:
            blobs.append(dot)

    print(f'Found {len(blobs)} dots.')
            
    return blobs


def grouping_dots(dot_list, dist=2):
    
    keep_dots = []
    
     # ensure dot_list is sorted by time
    frames = [one_dot.frame for one_dot in dot_list]
    sorted_i = np.argsort(frames) # get frame/time sorted indexes
    dot_list = np.array(dot_list)[sorted_i] # sort dot list in the same way -- earliest time to lastest time

    # {frame1: [list of dots], frame2: [list of dots]}
    dot_dict = {} # each frame is a dictonary key, value is list of dots found in that frame

    for ddot in dot_list: # for each dot observation found
        if ddot.frame not in dot_dict:
            dot_dict[ddot.frame] = [ddot] # if frame is not in dict, make new entry with dot observation
        else:
            dot_dict[ddot.frame].append(ddot) # if frame is already in dict, append dot observation
    
    frames_with_dots = np.sort(list(dot_dict.keys()))
    
    master_matched = [] # keep list of all matched dots -- ensure none are used twice
    
    for frame_curr in frames_with_dots:
    
        for one_dot_curr in dot_dict[frame_curr]: # for dot observations in the current frame
            
            if one_dot_curr in master_matched: # if dot has already been matched, skip
                continue

            ###
            def search_next_frame(curr_obs_list, frame_next): # function for checking if dot in next frame

                matches = {} # make dictionary of matched dots
                one_dot_curr = curr_obs_list[-1]

                for one_dot_next in dot_dict[frame_next]: # for dot observation in the next frame
                    
                    if one_dot_next not in master_matched: # if dot has not already been matched

                        is_same, dist = one_dot_curr.check_if_same_dot(one_dot_next)

                        if is_same: # if finds match at least once, append to match array
                            matches[dist] = one_dot_next
                            
                    else:
#                         print('already matched: ', one_dot_next)
                        continue

                if len(matches) > 0: # if found at least one match

                    # once finished looping through next frame, check which one has minimum distance
                    min_dist_obs = matches[np.nanmin(list(matches.keys()))]

                    # append minimum distance to observation list
                    curr_obs_list.append(min_dist_obs)

                    return True, curr_obs_list # return True, updated observation list

                else: # if no match found
                    return False, curr_obs_list # return False, unchanged observation list   
            ###

            searching = True
            observation_list = [one_dot_curr] # initialize observation list
            frame_next = frame_curr + 1 # initialize frame_next as the frame after frame_curr
            
#             print('trying to match: ', one_dot_curr)

            while searching:
                
                if frame_next in frames_with_dots:

                    found_match, observation_list = search_next_frame(observation_list, frame_next)

                    if found_match: # if match is found
                        frame_next += 1 # search next frame

                    else:
                        searching = False # if not match found, stop searching and break loop
                        
                else: # if frame_next does not contain any dots, stop searching
                    searching = False
    
            if len(observation_list) > 2: # if more than 2 observations of the dot

                for i in observation_list:
#                     print('   matched:', i)
                    master_matched.append(i) # for each matched observation, append to master_matched list

                keep_dots.append(dot(observation_list))
            
#             else:
#                 print('cannot match: ', observation_list[0])
        
    return keep_dots  



# def grouping_dots(dot_list, dist=2):
    
#     dot_grouped_list = []
    
#     # ensure dot_list is sorted by time
#     times = [one_dot.time.datetime.time() for one_dot in dot_list]
#     sorted_i = np.argsort(times) # get time sorted indexes
#     dot_list = np.array(dot_list)[sorted_i] # sort dot list in the same way -- earliest time to lastest time

#     matched_i = [] # list of matched observations indexes
#     for dot1_i, dot1 in enumerate(tqdm(dot_list)): # loop through dot observations

#         dot1_frame = dot1.frame # get frame of dot1
#         next_frame = dot1_frame + 1 # initialize next_frame variable starting with the frame after the first observation

#         if dot1_i not in matched_i: # if this dot has NOT already been matched

#             one_dot = dot([dot1]) # initialize dot object
#             matched_i.append(dot1_i) # append dot1 index to "matched_i" list

#             # matched_dots = [] # save matched dot observations
#             # matched_dists = []
#             # i_list = []

#             for dot2_i in np.arange(dot1_i+1, len(dot_list)): # loop through dot observations starting with the one after the dot1 observation

#                 dot2 = dot_list[dot2_i] # get dot2
#                 dot2_frame = dot2.frame # get frame of dot2

#                 if dot2_frame == next_frame: # if the observation to check is in the next frame

#                     if dot2_i not in matched_i: # double check this dot has not already been matched

#                         # check if match
#                         is_same_dot, diff_dist = dot1.check_if_same_dot(dot2, dist) # check if same dot, True if is

#                         # if match, append observation
#                         if is_same_dot: 

#                             # matched_dots.append(dot2)
#                             # matched_dists.append(diff_dist)
#                             # i_list.append(dot2_i)

#                             one_dot.add_obs(dot2) # add dot2 to one_dot object
#                             matched_i.append(dot2_i)

#                             #if the observation is a match, change "next_frame" to the frame after
#                             next_frame += 1
                            
#                 elif dot2_frame < next_frame: # if dot2_frame is before next_frame
#                     continue # continue to next observation

#                 # elif dot2_frame > next_frame: # observations sorted by time, so once observations in next frame won't have one before

#                     # if reached here, no longer have to keep searching! -- if obs matched, add it; if no obs matched, break loop

#                     # if len(matched_dots) > 0: # if there is at least one match

#                     #     min_i = np.argmin(matched_dists)

#                     #     print(f'dists: {matched_dists}')
#                     #     print(f'x: {[matched_dots[i].x for i in range(len(matched_dots))]}')
#                     #     print(f'y: {[matched_dots[i].y for i in range(len(matched_dots))]}')
#                     #     print(f'frames: {[matched_dots[i].frame for i in range(len(matched_dots))]}')

#                     #     one_dot.add_obs(matched_dots[min_i]) # add minimum distance obs to one_dot object
#                     #     matched_i.append(i_list[min_i]) # add obs index to matched_i list

#                     #     # reset minimum distance variables!
#                     #     matched_dots = [] # save matched dot observations
#                     #     matched_dists = []
#                     #     i_list = []


#                         # if the observation is a match, change "next_frame" to the frame after
#                         # next_frame += 1

#                     # else: # if there is no match in this frame, break the matching loop
#                     #     break # assume dot has disappeared and stop searching

#                 else:
#                     break
                    
#             if one_dot.get_lifetime() >= 10 * u.s: # if more than one observation of this dot -- EUI candence 5s so 10s lifetime ensures dot there for at least 2 frames, also 10s is IRIS cadence
#                 dot_grouped_list.append(one_dot)
#             else: # if not more than one obs of this dot, omit
#                 continue

#         else: # if this obs has already been matched, continue to next obs
#             continue

#     return dot_grouped_list

    
    # dots_checked_i = 0
    # while dots_checked_i < len(dot_list) - 1:
        
    #     dot1 = dot_list[dots_checked_i]
    #     one_dot = dot([ dot1 ]) # initialize dot object
        
    #     unmatched_frames = 0
    #     for dot2 in dot_list[dots_checked_i+1:]: # starting with the next entry in dot_list
            
    #         dots_checked_i += 1 # checked another dot so add to checked dot counter
        
    #         is_same_dot = dot1.check_if_same_dot(dot2, dist) # check if same dot, True if is
        
    #         if is_same_dot:
    #             one_dot.add_obs(dot2) # add dot2 to one_dot object

    #             # if match, reset unmatch counter to 0
    #             unmatched_frames = 0
                
    #         elif unmatched_frames > 2: # if two frames without a match with dot1, presume dot no longer exists
    #             # dots_checked_i =- 2 # consider the unmatched frames for the next dot
    #             break # break inner loop
                
    #         else:
    #             unmatched_frames += 1 # no match, add 1 to unmatched counter
    #             continue
        
    #     # if more than one observation of this dot -- EUI candence 5s so 10s lifetime ensures dot there for at least 2 frames, also 10s is IRIS cadence
    #     if one_dot.get_lifetime() > 10 * u.s:
    #         dot_grouped_list.append(one_dot)
    #     else:
    #         continue


def plot_all_dots(dot_list, figsize=(10,10), fontsize=12):
    """plotting all dot at their peak for all dots in the given dot list"""

    rand_i = np.arange(len(dot_list))

    # make rows of 10
    num_rows = int(np.ceil(len(dot_list) / 10))

    # setting up figure
    fig, axs = plt.subplots(figsize=figsize)
    axs.axis('off')

    for plotting_i, i in enumerate(tqdm(rand_i)):
        one_dot = dot_list[i]

        try:
            one_dot.fit_gaussians()
            peak_index = one_dot.get_peak_i()

            one_image = one_dot.dot_obs_list[peak_index].image
            one_data = np.array(one_image.data)

            region_r = 4
            blob_x = one_dot.get_xs()[peak_index]
            blob_y = one_dot.get_ys()[peak_index]

            # blob region
            x1, x2, y1, y2 = int(blob_x-region_r), int(blob_x+region_r), int(blob_y-region_r), int(blob_y+region_r)
            one_data = one_data[y1:y2, x1:x2]

            # plotting
            ax = fig.add_subplot(num_rows, 10, 1+plotting_i, projection=one_image.wcs)
            ax.axis('off')
            ax.set_title(f'Dot {i}', fontsize=fontsize)
            ax.imshow(one_data, cmap=one_image.cmap)
            # ax.set_xlim(x1, x2)
            # ax.set_ylim(y1, y2)

        except:
            continue
        
    plt.tight_layout()
    plt.show()

    return fig


def check_if_same_dot(dot1, dot2):

    xs1 = dot1.get_xs() # if present in the same frames and at the same x positions, assume the same dot
    frames1 = dot1.get_frames()

    xs2 = dot2.get_xs()
    frames2 = dot2.get_frames()
    
    if (list(xs1) == list(xs2) and frames1 == frames2) or (xs2.all() in xs1 and np.array(frames2).all() in np.array(frames1)) or (xs1.all() in xs2 and np.array(frames1).all() in np.array(frames2)):
        return True
    else:
        return False



def plot_all_dots_LC(eui_dots, iris_dots, figsize=(10,10), fontsize=12, time_thresh=20, colors = ['#FFC107', '#D81B60']): 
    """ time_thresh: (seconds) time between iris and eui intensity peaks to consider them at the same time
    """
    rand_i = np.arange(len(eui_dots))

    # make rows of 10
    num_rows = int(np.ceil(len(eui_dots) / 10))

    # setting up figure
    fig, axs = plt.subplots(figsize=figsize)
    axs.axis('off')
    
    same_time_peak = 0

    # FIND THE IRIS DOTS THAT ARE MATCHED TO MORE THAN ONE EUI DOT
    iris_dots_repeat_i = []
    for i, one_dot in enumerate(iris_dots):
        
        # one_xs = one_dot.get_xs()
        # one_frames = one_dot.get_frames()
        
        for check_i, check_dot in enumerate(iris_dots[i+1:]):
            
            true_i = i+1 + check_i

            is_same = check_if_same_dot(one_dot, check_dot)

            if is_same:
                iris_dots_repeat_i.append(i)
                iris_dots_repeat_i.append(true_i)

            # check_xs = check_dot.get_xs()
            # check_frames = check_dot.get_frames()
            
            # if (list(check_xs) == list(one_xs) and check_frames == one_frames) or (check_xs.all() in one_xs and np.array(check_frames).all() in np.array(one_frames)) or (one_xs.all() in check_xs and np.array(one_frames).all() in np.array(check_frames)):
            #     iris_dots_repeat_i.append(i)
            #     iris_dots_repeat_i.append(true_i)


    for plotting_i, dot_i in enumerate(tqdm(rand_i)):
        
        ax = fig.add_subplot(num_rows, 10, 1+plotting_i)
            
        ax.set_title(f'Dot {dot_i}', fontsize=fontsize)
        ax.tick_params(axis='x', labelrotation = -45)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)
            
        eui_times = eui_dots[dot_i].get_earth_times()
        eui_intensity = eui_dots[dot_i].get_enhancement()
        norm_eui_intensity = eui_intensity / np.nanmax(eui_intensity)
        ax.plot(eui_times, norm_eui_intensity, c=colors[0])
        ax.scatter(eui_times, norm_eui_intensity, label='EUI', c=colors[0], s=10)

        iris_times = iris_dots[dot_i].get_earth_times()
        iris_intensity = iris_dots[dot_i].get_enhancement()
        norm_iris_intensity = iris_intensity / np.nanmax(iris_intensity)
        if dot_i in iris_dots_repeat_i:
            ax.plot(iris_times, norm_iris_intensity, c=colors[1], linestyle=':', alpha=0.5) # if repeat iris dot, make dotted line
        else:
            ax.plot(iris_times, norm_iris_intensity, c=colors[1], alpha=0.5)
        ax.scatter(iris_times, norm_iris_intensity, label='IRIS', c=colors[1], s=10)
        
        # get when iris at peak intensity
        peak_iris = iris_dots[dot_i].get_peak_i()
        peak_eui = eui_dots[dot_i].get_peak_i()
        if not np.isnan(peak_iris): # if it has a corresponding iris dot
            iris_peak_time = iris_times[peak_iris]
            eui_peak_time = eui_times[peak_eui]
            if abs(iris_peak_time - eui_peak_time) / np.timedelta64(1, 's') <= time_thresh:
                ax.set_title(f'Dot {dot_i}', fontsize=fontsize, fontweight='bold')
                same_time_peak += 1

        
        min_start = np.min([list(eui_times) + list(iris_times)])
        max_end = np.max([list(eui_times) + list(iris_times)])
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        ax.set_ylim(0, 1.1)
        
    print(f'IRIS + EUI dots peaking within 20s: {same_time_peak} dots, {same_time_peak/(len(iris_dots)) * 100:.2f}%')
    print(f'There are {len(iris_dots_repeat_i)} EUI dots that use the same IRIS dots')
        
    plt.tight_layout()
    plt.show()
    
    return fig


def convert_pix_to_dist(image, pix_val):
    pixscale_x = image.scale[0] # deg/pix in x
    pixscale_y = image.scale[1] # deg/pix in y
    d_sun = image.dsun # Observer distance from the center of the Sun
    return (((pix_val*u.pix * pixscale_x).to(u.rad) * d_sun) / (u.rad)).to(u.km).value # km
        
