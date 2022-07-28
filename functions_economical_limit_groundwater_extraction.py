#
# Global Economic Limits of Groundwater for Irrigation
#
# Copyright (c) L.P.H. (Rens) van Beek / Marc F.P. Bierkens 2018-2022
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###########
# modules #
###########

import numpy as np
import pcraster as pcr

#############
# functions #
#############

# Ancillary functions to read the input from tables.
# Eventually, this may be replaced by more generic functions and ones that
# explicitly consider time.

def pcr_get_map_value( \
        pcrfield, 
        location = 1, \
        pcrfunc = 'areaaverage', \
        mv = -999.9, \
        format_str = '%.3f', \
        message_str = '', \
        test_verbose = False):

    '''
pcr_print_cell_value: generic function to return the value from a PCRaster field \
at a specified location. The location specified may be a cell number (default), \
a tuple containing the row, column number or a classified PCRaster map identifying \
points or locations. Output in that case is then dependent on the PCRaster area ... \
function specified and consists of a dictionary of key, value pairs; \
otherwise, the value is returned.

    Input:
    ======
    pcrfield:           PCRaster field for which the values are extracted; cast
                        as scalar in the function;
    location:           optional argument, default value is 1 (top-left cell);
                        the location specified may be a cell number (default),
                        or a list of locations or tuples containing the row,
                        column number or a classified PCRaster map identifying
                        points or areas;
    pcrfunc:            PCRaster function or string of the function name, \
                        either being the areaaverage, areatotal, areaminimum or
                        areamaximum. Default value is areaaverage; if the loc-
                        ation is not a map, this argument is ignored;
    mv:                 standard missing value identifier, default = -999.9;
    format_str:         string of format, default is %.3f;
    message_str:        text that is added to the returned message string,
                        default value is an empty string;
    test_verbose:       boolean variable, that will print the message string to
                        the screen if True; default value is False.

    Output:
    =======
    value:              value (float) or dictionary of floats retrieved from
                        the map on the basis of the specified location or loc-
                        ations. A single value is returned if only a single
                        location is included;
    message_str:        the composed message_str.

'''

    # initialize
    # set the value
    value = {1: None}
    
    # cast all variables as fields
    pcrfield = pcr.spatial(pcr.scalar(pcrfield))

    # test if the location is a PCRaster field
    if type(location) is pcr.Field:

        # set the function
        if type(pcrfunc) is str and len(pcrfunc) > 0:
            pcrfunc = getattr(pcr, pcrfunc)
        else:
            pcrfunc = pcr.areaaverage
            
        pcrfield = pcr.cover(pcrfunc(pcr.spatial(pcrfield), \
            location), mv)

        # retrieve the locations
        loc_a = pcr.pcr2numpy(location, 0)
        loc_ids = np.unique(loc_a)
        loc_ids= (loc_ids[loc_ids != 0]).tolist()
        loc_ids.sort()
        
        # insert the values
        val_a = pcr.pcr2numpy(pcrfield, 0)
        for loc_id in loc_ids:
            value[loc_id] = val_a[loc_a == loc_id][0]
            message_str = str.join('\n', (message_str, \
                '%3d: %s' % (loc_id, \
                     format_str % value[loc_id])))

    # location is not a PCRaster field
    else:
        
        # cast location as a list
        if not type(location) is list:
            location = [location]

        # extract the location iteratively
        for loc_id in range(len(location)):
            if type(location[loc_id]) is tuple:
                # get the value
                valx, valid = pcr.cellvalue(pcrfield, \
                        location[loc_id][0], location[loc_id][1])
            else:
                # get the value
                valx, valid = pcr.cellvalue(pcrfield, \
                        location[loc_id])
            # value returned, process
            if not valid:
                valx = mv
            # add valx to value and update the message string
            value[loc_id + 1] = valx
            message_str = str.join('\n', (message_str, \
                '%3d: %s' % (loc_id + 1, format_str % value[loc_id + 1])))

    # all data added, print the message string if test_verbose is True
    if test_verbose:
        print (message_str)

    # return the direct value
    if len(value) == 1:
        key = list(value.keys())[0]
        value = value[key]

    # return the value and message_str
    return value, message_str

# read_water_productivity ######################################################
def read_water_productivity(datafile, selectedcrops, testverbose= False):
    '''
Reads in the water use efficiency from a file in which \
the water use efficiency in kg per m3 of evapotranspirated water is \
stored per crop type.

    Input:
    ======
    datafile:               name of data file, including full path
                            comments can be stored when preceded by '#'
                            data are stored with one entry per crop per column on one row
    selectedcrops:          list of selected crop ids, corresponding with the row number
    testverbose:            optional; prints selected output

    Output:
    =======
    water_ productivity:    water productivity for each crop [kg*m^-3]

'''
    
    #-start
    print (str.join(' ', ('\t*** NOTE: this may be replaced by more generic function', \
        'that explicitly considers time.')))
    if testverbose:
        print ('\tinitializing water productivity per crop')

    water_productivity = {}
    f= open(datafile)
    try:
        for line in f:
            if line[0]!= '#':
                rawlist= line.split(";")
                for icnt in range(len(rawlist)):
                    crop_id = icnt + 1
                    if crop_id in selectedcrops:
                        water_productivity[crop_id]=  float(rawlist[icnt])              
    finally:
        f.close()

    # patch any missing values
    for crop_id in selectedcrops:
        if crop_id not in water_productivity.keys():
            water_productivity[crop_id] = 0.00
    
    if testverbose:
        print ('\twater use efficiency read from file %s' % datafile)
        for crop_id in selectedcrops:
            print ('\t%16s: %.3f' % (crop_id, \
                water_productivity[crop_id]))

    # return the water productivity
    return water_productivity

# read_dictionary_from_table ###################################################
def read_dictionary_from_table(datafile, key_column, data_column, \
        testverbose = False):
    '''
Reads in miscellaneous information per key (e.g., country) from a file with \
tabulated data where the keys are stored across the rows and the information is \
contained by a single column.

    Input:
    ======
    datafile:               name of data file, including full path
                            comments can be stored when preceded by '#' and
                            information is expected to be organized with the
                            variables across the columns and the identifiers (keys)
                            along the rows;
    key_column:             column in data table to use as key in the resulting dictionary
    data_column:            column in data table to use as index to retrieve info
                            for he resulting dictionary
    testverbose:            optional; prints selected output

    Output:
    =====
    data_info:               information from the data table stored as a dictionary

'''

    #-start
    print (str.join(' ', ('\t*** NOTE: this may be replaced by more generic function', \
        'that explicitly considers time.')))
    if testverbose:
        print ('\tretrieving information from table %s' % datafile)

    data_info= {}
    f = open(datafile)
    try:
        for line in f:
            if line[0]!= '#':
                rawlist = line.split(";")
                # get the key
                try:
                    key = int(rawlist[key_column])
                except:
                    key = rawlist[key_column].strip()
                # add the value
                try:
                    data_info[key]= float(rawlist[data_column].strip())
                except:
                    data_info[key]= rawlist[data_column].strip()
    finally:
        f.close()
    if testverbose:
        print ('\t%d entries were retrieved from data column %d' %\
            (len(data_info), data_column))
        for key in data_info.keys():
            # check and act on type
            if type(data_info[key]) is float:
                m_str = '%.1f' % data_info[key]
            else:
                m_str = str(data_info[key])
            # print the information
            print ('\t%16s: %s' % (str(key), \
                data_info[key]))

    # return the dictionary
    return data_info

# map_info_to_pcr ##############################################################
def map_info_to_pcr(data_info, key_map, pcr_data_type):

    '''
Reads in miscellaneous information per key (e.g., country) from a dictionary \
where the keys are stored across the rows and the information is contained by \
a single column and returns the corresponding map.

    Input:
    ======
    data_info:              dictonary that contains key, value pairs that are
                            mapped to the key map specified;
    key_map:                a classified PCRaster field which entries correspond
                            with those of the keys;
    pcr_data_type:          PCRaster function to cast data type (e.g., pcr.scalar).
    output:
    =====
    result_map:             information from the data table stored as a scalar
                            PCRaster field.

'''

    #-insert values in map
    result_map = pcr_data_type(0)
    
    # iterate over the information and update the resulting map
    for key, value in data_info.items():
        result_map = pcr.ifthenelse(key_map == key, pcr_data_type(value), \
                result_map)
    
    # all values inserted
    # clip to mask
    result_map = pcr.ifthen(pcr.defined(key_map), result_map)

    # return the resulting map
    return result_map

# map_table_info_to_pcr ########################################################
def map_table_info_to_pcr(datafile, key_map, key_column, data_column, \
    pcr_data_type, testverbose= False):

    '''
Reads in miscellaneous information per key (e.g., country) from a file  where\
the keys are stored across the rows and the information is contained by \
a single column and returns the corresponding map.

    Input:
    ======
    datafile:               name of data file, including full path
                            comments can be stored when preceded by '#' and
                            information is expected to be organized with the
                            variables across the columns and the identifiers (keys)
                            along the rows;
    key_map:                a classified PCRaster field which entries correspond
                            with those of the key column;
    key_column:             column in data table to use as key to link data to
                            the classified map;
    data_column:            column in data table to use as index to look up data
                            in the table and link them to the keys;
    pcr_data_type:          PCRaster function to cast data type (e.g., pcr.scalar);
    testverbose:            optional; prints selected output

    output:
    =====
    result_map:             information from the data table stored as a scalar
                            PCRaster field.

'''
    
    if testverbose:
        print ('\tinitializing map from table %s for data column %d using column %d as key' % \
            (datafile, data_column, key_column))

    #-get data info from table
    data_info = read_dictionary_from_table(datafile = datafile, \
        key_column = key_column, data_column = data_column, \
        testverbose = testverbose)

    #-insert values in map
    result_map = map_info_to_pcr(data_info, key_map, pcr_data_type)

    if testverbose:
        print ('\tmap from file initialized')
    
    # return the resulting map
    return result_map

# sum_list #####################################################################

def sum_list(list_in):
    '''
sum_list: generic function that returns the sum of a list of which \
the entries can be summed (integers, floats, arrays etc.).

'''
    #-sums all entries in a list-like type
    return sum(list_in)

# def multiply_dicts ###########################################################

def sum_dicts(a_dict, b_dict, halt_on_key_error = False):
    '''
sum_dicts: generic function that returns a dictionary that uses the inter-\
section of the keys of dictionaries a and b and the sum of their entries as \
values. This assumes that the  entries can be summed directly (integers, \
floats, arrays etc.) and missing values are not explicitly handled by the \
function.

The function can produce a warning or halt if keys of the dictionaries a and b \
do not match if halt_on_key_error is set respecively False or True \
(default: False).

'''

    # create an empty dictionary
    c_dict = {}
    
    # check on keys
    if len(a_dict) != len(b_dict):
        message_str = 'length of keys of dictionaries (%d, %d) do not match' % \
            (len(a_dict), len(b_dict))
        if halt_on_key_error:
            sys.exit('error: %s' % message_str)
        else:
            print ('warning: %s' % message_str)
    
    # extract common keys
    common_keys = []
    if len(a_dict) >= len(b_dict):
        keys = list(a_dict.keys())
    else:
        keys = list(b_dict.keys())
    for key in keys:
        if key in a_dict.keys() and key in b_dict.keys():
            common_keys.append(key)
    
    # sort the common keys
    common_keys.sort()
    
    # get the values
    for key in common_keys:
        c_dict[key] = a_dict[key] + b_dict[key]

    # and return the resulting dictionary
    return c_dict 

# def multiply_dicts ###########################################################

def multiply_dicts(a_dict, b_dict, halt_on_key_error = False):
    '''
multiply_dicts: generic function that returns a dictionary that uses the inter-\
section of the keys of dictionaries a and b and the product of their entries as \
values. This assumes that the  entries can be multiplied directly (integers, \
floats, arrays etc.) and missing values are not explicitly handled by the \
function.

The function can produce a warning or halt if keys of the dictionaries a and b \
do not match if halt_on_key_error is set respecively False or True \
(default: False).

'''

    # create an empty dictionary
    c_dict = {}
    
    # check on keys
    if len(a_dict) != len(b_dict):
        message_str = 'length of keys of dictionaries (%d, %d) do not match' % \
            (len(a_dict), len(b_dict))
        if halt_on_key_error:
            sys.exit('error: %s' % message_str)
        else:
            print ('warning: %s' % message_str)
    
    # extract common keys
    common_keys = []
    if len(a_dict) >= len(b_dict):
        keys = list(a_dict.keys())
    else:
        keys = list(b_dict.keys())
    for key in keys:
        if key in a_dict.keys() and key in b_dict.keys():
            common_keys.append(key)
    
    # sort the common keys
    common_keys.sort()
    
    # get the values
    for key in common_keys:
        c_dict[key] = a_dict[key] * b_dict[key]

    # and return the resulting dictionary
    return c_dict 

# redistribute_pcr_field_by_class_mean #########################################

def redistribute_pcr_field_by_class_mean( \
    value, class_weights, class_values, area_id):

    '''
redistribute_population_by_cultivated_area: function that redistributes a PCRaster
value field on the basis of the area-averaged weight into totals for the area and
then redistributes them on the basis of the class values and the area-averaged ratio.

'''

    # set the initial value
    value_ini = value
    
    # start by computing the weights
    total_weight = sum(class_weights)
    
    # and the total value
    total_value_ini = pcr.areatotal(value_ini, area_id)
    
    # set the area averaged mean
    area_means = [None for class_weight in class_weights]
    values = [None for class_weight in class_weights]

    # iterate over the class weights
    for icnt in range(len(class_weights)):
        
        # set the value for standardization
        class_weight = class_weights[icnt]
        mask = pcr.defined(class_weight) & (total_weight > 0)
        
        # class_weight: this is a fraction per area, summing to unity
        class_weight = pcr.areatotal(class_weight, area_id) / pcr.areatotal(total_weight,area_id)
        class_weights[icnt] = class_weight
        
        # area-averaged ratio between the class values and the fraction of the input
        # value based on the weights; start by resetting the mask
        mask = class_weight * total_value_ini > 0
        area_means[icnt] = pcr.ifthenelse(mask, \
            pcr.areatotal(class_values[icnt], area_id) / \
            (class_weight * total_value_ini), 0)

        # and set the values
        values[icnt] = pcr.ifthen(pcr.defined(value_ini) & \
            pcr.defined(class_values[icnt]), \
                pcr.cover(class_values[icnt] / area_means[icnt], 0))

    # get value and total value
    value = sum(values)
    total_value = pcr.areatotal(value, area_id)
    
    # reset the areas that became zero
    value = pcr.ifthenelse(total_value > 0, value, value_ini)
    total_value = pcr.areatotal(value, area_id)

    # return the final value
    return value

# compute_groundwater_depth_from_extraction ####################################

def compute_groundwaterdepth_from_extraction(groundwaterextraction, \
        groundwaterdepth, cellarea, \
        specificyields, layerdepths, mv_id = -999.9, testverbose = False):
    
    # set the initial top of the layer
    layertop = pcr.spatial(pcr.scalar(0))
    
    # set a test location
    if testverbose:
        test_depth = 243.0
        test_location = pcr.cover(pcr.order(pcr.abs(layerdepths[0] - \
                                  pcr.ifthen(groundwaterextraction > 0, \
                                             pcr.scalar(test_depth)))), 0) == 1
    
    # iterate over the layers and extract water
    for ix in range(len(specificyields)):
        
        # get the specific yield and bottom
        specificyield = specificyields[ix]
        layerbottom = layerdepths[ix]
        # replace the bottom layer depth if it is the latest value
        # and the value is a missing value
        if ix == len(specificyields) - 1:
            
            # remove the depth
            layerbottom = pcr.ifthenelse(layerbottom == mv_id, \
                                         2 * (groundwaterdepth + groundwaterextraction /\
                                         (specificyield * cellarea)), layerbottom)
        
        # get the effective depth
        effectivedepth = pcr.max(0, layerbottom - groundwaterdepth)

        if testverbose == True:
            print ('\n* Processing layer %d' % ix)
            print (pcr_get_map_value(cellarea, test_location, \
                   message_str = 'cell area:')[1])           
            print (pcr_get_map_value(layerbottom, test_location, \
                   message_str = 'layer depth:')[1])
            print (pcr_get_map_value(groundwaterdepth, test_location, \
                   message_str = 'groundwater depth:')[1])
            print (pcr_get_map_value(specificyield, test_location, \
                   message_str = 'specific yield:')[1])
            print (pcr_get_map_value(effectivedepth, test_location, \
                   message_str = 'effective depth:')[1])

        # extractable water is the minimum of what is available in the layer and
        # what is needed
        extractablewater = effectivedepth * specificyield * cellarea
        extractablewater = pcr.min(extractablewater,  groundwaterextraction)

        if testverbose == True:
            print (pcr_get_map_value(groundwaterextraction, test_location, \
                   message_str = 'groundwater extraction:')[1])
            print (pcr_get_map_value(extractablewater, test_location, \
                   message_str = 'available abstraction:')[1])


        # and update the groundwater depth and extraction
        groundwaterdepth = groundwaterdepth + extractablewater / \
                           (specificyield * cellarea)
        groundwaterextraction = pcr.max(0, groundwaterextraction - extractablewater)
        
        if testverbose == True:
            print (pcr_get_map_value(groundwaterdepth, test_location, \
                   message_str = 'groundwater depth:')[1])                  
            print (pcr_get_map_value(groundwaterextraction, test_location, \
                   message_str = 'groundwater extraction:')[1])       
        
        # reset the top of the layer
        layertop = layerbottom

    # return the groundwater depth
    return groundwaterdepth

# compute_overpressure_sprinkers ###############################################
def compute_overpressure_sprinklers(irr_radius, \
    grav_acc = 9.81, rho_water = 1000.0, cor_offset = 0.00, cor_factor = 1.00):
    '''compute_overpressure_sprinklers: function to compute the overpressure in [Pa] for \
sprinkler irrigation systems under the assumption that the water is pumped through \
the system without resistance and does not experience air drag on its path through \
the air. Furthermore, the jet is assumed to be ejected under an angle of 45 degrees \
with the horizontal at an elevation of zero metres on an otherwise plane field and cover \
a circular irrigated area of radius R.

    Input:
    ======
    irr_radius:             radius [m] of the circle that is irrigated by the jet
                            that is ejected by the sprinkler under the applied
                            overpressure;
    grav_acc:               gravitational acceleration, set to 9.81 [m/s2]
                            by default
    rho_water:              density of water, set to 1000 [kg/m3] by default
    cor_offset:             additive correction factor, default 0.00 [-];
    cor_factor:             multiplicative correction factor, default 1.00 [-],
                            correction factors are only applied to the over-
                            pressure, not to the velocity.
    
    Output:
    =======
    overpressure:           overpressure in [Pa] to create a velocity of the jet
                            that is large enough to cover the distance of the
                            irrigated circle of radius, R, when ejected under an
                            angle of 45 degrees;
    velocity:               velocity of ejection corresponding with the over-
                            pressure in [m/s].
   
'''
    # compute the overpressure and the velocity
    velocity = (grav_acc * irr_radius) ** 0.5
    overpressure=  0.25 * velocity ** 2 * rho_water
    overpressure = cor_offset + cor_factor * overpressure
    
    # return the values
    return overpressure, velocity

def compute_pumping_efficiency(groundwaterdepth, initial_pumping_efficiency = 0.95, \
        final_pumping_efficiency = 0.60, depth_norm = 180.0):
    
    '''dummy: returns the efficiency that approximates 0.8 at 100 m and ends at 0.6 at 800 m.'''
    
    initial_pumping_efficiency = pcr.min(1.0, pcr.scalar(initial_pumping_efficiency))

    pumping_efficiency = initial_pumping_efficiency + \
        (final_pumping_efficiency - initial_pumping_efficiency) * \
        (1.00 - pcr.exp(-groundwaterdepth / depth_norm))
    
    return pumping_efficiency
    

#~ # compute_deep_well_turbine_efficiency #########################################
    #~ 
#~ def compute_deep_well_turbine_efficiency(maximum_efficiency = 0.80, \
        #~ optimum_):
    #~ '''
#~ compute_deep_well_turbine_efficiency: function that returns the efficiency.
#~ 
#~ '''
    
    

#~ # readMappedInfoFromTable ######################################################
#~ def readMappedInfoFromTable(dataFile, keyMap, keyColumn, dataColumn, pcrDataFunc, verbose= False):
    #~ 
    #~ '''
#~ Reads in miscellaneous information per key (e.g., country) from a file  where\
#~ the keys are stored across the rows and the information is contained by \
#~ a single column and returns the corresponding map.
#~ 
    #~ Input:
    #~ ======
    #~ dataFile:               name of data file, including full path
                            #~ comments can be stored when preceded by '#'
                            #~ data are stored with one entry per crop per column on one row
    #~ keyMap:                 a PCRaster map which entries correspond with those of the keyed column
    #~ keyColumn:              column in data table to use as key in the resulting dictionary
    #~ dataColumn:             column in data table to use as index to retrieve info
                            #~ for he resulting dictionary
    #~ pcrDataFunc:            PCRaster function to cast data type
    #~ verbose:                optional; prints selected output
#~ 
    #~ Output:
    #~ =====
    #~ mapInfo:               information from the data table stored as a PCRaster map
#~ 
#~ '''
    #~ print '\t*** TEMPORARY: replace with generic function to read data ***'
    #~ if verbose:
        #~ print '\tinitializing map from table info stored as file'
#~ 
    #~ #-get data info from table
    #~ dataInfo= readDictInfoFromTable(dataFile, keyColumn, dataColumn, verbose= verbose)
    #~ 
    #~ #-insert values in map
    #~ mapInfo= pcrDataFunc(0)
    #~ for key, value in dataInfo.iteritems():
        #~ mapInfo= pcr.ifthenelse(keyMap == key, pcrDataFunc(value), mapInfo)
    #~ # all values inserted
    #~ # clip to mask
    #~ mapInfo= pcr.ifthen(pcr.defined(keyMap), mapInfo)
#~ 
    #~ if verbose:
        #~ print '\tmap from file initialized'
    #~ 
    #~ return mapInfo
    #~ 

# compute energy requirement with depth
def compute_total_energy_well_extraction(total_pumped_water_volume, \
        potential_head_difference, pumping_efficiency,\
        rhoWater= 1004.0, gravAcc= 9.81):

    '''
Computes the total energy required in MJ to pump the required water volume \
for the specified static head and the corresponding pumping efficiency. \
Takes also preset variables for the density of \
water, and the gravitational acceleration constant.

    Input:
    ======
    total_pumped_water_volume:    total pumped water requirement, m3/year
    potential_head_difference:    head difference in potenial energy field, m
    pumping_efficiency:           dimensionless, ratio of efficiently pumped
                                  water volume (< 1, typically 0.6-0.8)
    rhoWater:                     optional: density of water, default 1004 kg/m3
    gravAcc:                      optional gravitational acceleration,
                                  default 9.81 m/s2

    Output:
    =======
    total_energy:                 total energy for volume water pumped for the
                                  current well depth in MJ

''' 
    
    return 1.0e-6 * (total_pumped_water_volume * rhoWater * gravAcc * \
            potential_head_difference) / pumping_efficiency


# compute_specific_energy_cost_with_depth ######################################
def compute_specific_energy_cost_withdrawal(total_pumped_water_volume, \
        potential_head_difference, electricity_price, pumping_efficiency,\
        rhoWater= 1004.0, gravAcc= 9.81, conversion_factor= 3.6e6):
    '''
Computes the total energy costs in kWh to pump the required water volume \
for the specified static head and the corresponding pumping efficiency \
and the electricity price.
Takes also preset variables for the density of \
water, the gravitational acceleration constant and a conversion factor.

    Input:
    ======
    total_pumped_water_volume:    total pumped water requirement, m3/year
    potential_head_difference:    head difference in potenial energy field, m
    electricity_price:            price of electricity in $ per kWh
    pumping_efficiency:           dimensionless, ratio of efficiently pumped
                                  water volume (< 1, typically 0.6-0.8)
    rhoWater:                     optional: density of water, default 1004 kg/m3
    gravAcc:                      optional gravitational acceleration,
                                  default 9.81 m/s2
    conversion_factor:             optional, conversion factor to obtain correct
                                  units ($ per m depth of the volume pumped),
                                  default 3.6e6

    Output:
    =======
    specificEnergyCost:           energy cost for volume water pumped
                                  in $ per well.

''' 

    energy_requirement = compute_total_energy_well_extraction(total_pumped_water_volume, \
        potential_head_difference, pumping_efficiency,\
        rhoWater, gravAcc)

    return electricity_price * energy_requirement * (1.0e6 / conversion_factor)

# compute_well_costs_with_depth ################################################
def compute_well_costs_with_depth(construction_cost_per_m, labour_cost_per_m, \
        well_depth, number_wells):
    '''
''' 

    return number_wells * well_depth * \
        (construction_cost_per_m + labour_cost_per_m)

################################################################################
    #~ 
#~ # computeMaximumProfitableWellDepth ############################################
#~ def computeMaximumProfitableWellDepth(totalIncome, specificLabourCost,\
    #~ specificMaterialCost, specificEnergyCost, wellWaterDepth):
    #~ '''
#~ 
#~ Computes the maximum profitable well depth, the depth where the total annual \
#~ costs equal the annual income from the well.
#~ Investment costs apply to the construction of the well to its full depth, \
#~ energy costs apply to the depth of the static head inside the well, \
#~ the difference has to be specified.
#~ 
    #~ Input:
    #~ ======
    #~ totalIncome:              total income per year in $       
    #~ specificLabourCost:       specific labour cost in $ per metre
    #~ specificMaterialCost:     specific material cost in $ per metre
                              #~ labour and material costs are investments
    #~ specificEnergyCost:       specific energy cost in $ per metre
    #~ wellWaterDepth:           depth of water in the well, metre
    #~ 
    #~ Output:
    #~ =======
    #~ totalWellDepth:           total depth to the bottom of the well, in m
    #~ wellWaterdepth:           wellWaterDepth, height of water in well, in m
    #~ conditionMask:            nominal mask, specifying conditions:
                              #~ 0: no irrigation, 1: economic limit
                              #~ of costs, 2: economic limit of income
    #~ 
#~ '''
#~ 
    #~ #-initialize mask of conditions: 0: no irrigation, 1: economic limit
    #~ # of costs, 2: economic limit of income 
    #~ conditionMask= pcr.nominal(0)
#~ 
    #~ #-start by computing the costs that are in excess of the well water depth
    #~ offSetCost= wellWaterDepth*(specificLabourCost+specificMaterialCost)
    #~ conditionMask= pcr.ifthenelse(totalIncome > 0,\
        #~ pcr.ifthenelse(offSetCost > totalIncome, pcr.nominal(1), pcr.nominal(2)),\
        #~ conditionMask)
    #~ offSetCost= pcr.min(offSetCost, totalIncome)
    #~ wellWaterDepth= offSetCost/(specificLabourCost+specificMaterialCost)
    #~ #-get remainder of well depth below the water table
    #~ totalWellDepth= wellWaterDepth+pcr.max(0, totalIncome-offSetCost)/(\
        #~ specificEnergyCost+specificLabourCost+specificMaterialCost)
    #~ 
    #~ #-return the total well depth, well water depth and condition mask
    #~ return totalWellDepth, wellWaterDepth, conditionMask
#~ 
#~ # redistributePopulation #######################################################
#~ def redistributePopulation(distributedValue, assignedValue, areaID):
#~ 
    #~ #-redistributes the scalar PCRaster map
    #~ # distributeValue over the selected area identified
    #~ # by the areaID to the assignValue; the surplus is determined
    #~ # and redistributed by the standardized value of assignValue
    #~ # on the basis of the area-averaged density
#~ 
    #~ #-copy the distributedValue
    #~ newDistributedValue= distributedValue
    #~ #-get masks where distributedValue and assignedValue
    #~ # are matched and where not
    #~ matchMask= pcr.cover((newDistributedValue > 0) & (assignedValue > 0), 0)
    #~ mismatchMask= pcr.cover((newDistributedValue == 0) & (assignedValue > 0), 0)
    #~ #-get target density
    #~ targetDensity= pcr.cover(pcr.areatotal(pcr.ifthen(matchMask, newDistributedValue), areaID)/\
    #~ pcr.areatotal(pcr.ifthen(matchMask, assignedValue), areaID), 0)
    #~ #-get surplus per cell and total, subtract the surplus from the cell value
    #~ # also, get the weight, so any remaining surplus can be reassigned
    #~ surplus= pcr.max(0, newDistributedValue-pcr.ifthenelse(matchMask, assignedValue, 0)*targetDensity)
    #~ surplusTotal= pcr.areatotal(surplus, areaID)
    #~ newDistributedValue= pcr.max(0, newDistributedValue-surplus)
    #~ surplusWeight= pcr.cover(surplus/surplusTotal, 0)
    #~ #-for the areas that have a mismatch, find the new value to add and the corresponding total
    #~ reassignWeight= pcr.cover(pcr.ifthenelse(mismatchMask, assignedValue, 0)/\
        #~ pcr.areatotal(pcr.ifthenelse(mismatchMask, assignedValue, 0), areaID), 0)
    #~ addition= pcr.cover(pcr.min(reassignWeight*surplusTotal, targetDensity*pcr.ifthenelse(mismatchMask, assignedValue, 0)), 0)
    #~ additionTotal= pcr.areatotal(addition, areaID)
    #~ #-get remainder
    #~ remainder= surplusWeight*pcr.max(0, surplusTotal-additionTotal)
    #~ #-update values and get mask of new condition
    #~ newDistributedValue+= addition+remainder
    #~ conditionMask= (newDistributedValue >= (assignedValue*targetDensity)) & (assignedValue > 0)
    #~ #-return the values
    #~ return newDistributedValue, conditionMask
#~ 
#~ # end of functions #############################################################
#~ 

#=========================================================================
#~ # compute the net income per capita per country
#~ country_population = pcr.areatotal(agricultural_population, countries)
#~ net_costs_per_capita  = pcr.areatotal(total_food_costs, countries) /\
    #~ country_population
#~ net_income_per_capita = pcr.areatotal(total_income - total_food_costs,\
    #~ countries) / country_population
#~ 
#~ # get the carrying capacity
#~ # no update will be performed if the net income per capita is negative
#~ target_population = pcr.ifthenelse(net_income_per_capita > 0, \
    #~ total_income / (net_income_per_capita + net_costs_per_capita), \
    #~ agricultural_population)
#~ 
#~ # get the total per country of the people that can be accommodated and
#~ # that is in surplus
#~ population_surplus = pcr.areatotal(pcr.max(0.0, \
    #~ agricultural_population - target_population), countries)
#~ loss = pcr.max(0.0, \
    #~ agricultural_population - target_population)
#~ gain = pcr.max(0.0, \
    #~ target_population - agricultural_population)
#~ agricultural_population = pcr.max(0.0, agricultural_population + \
    #~ gain - loss)
#~ agricultural_population = agricultural_population * \
    #~ pcr.cover(1.0, country_population / \
    #~ pcr.areatotal(agricultural_population, countries))
#~ 
#~ # check on convergence
#~ difference_map = pcr.areatotal(\
    #~ pcr.abs(agricultural_population - target_population), \
    #~ countries) / country_population
#~ 
#~ if pcr.cellvalue(pcr.mapmaximum(difference_map), 1)[0] < 1.0e-6:
                #~ update_population = False
#===============================================================================
