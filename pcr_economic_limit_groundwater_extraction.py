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
import os
import sys
import datetime
import logging

import numpy as np
import pcraster as pcr
import pcraster.framework as pcrm

sys.path.insert(0,'/home/rens/Scripts/pylib')

from spatialDataSet2PCR import spatialAttributes, setClone, spatialDataSet
from ncRecipes_fixed import getNCDates, writeField
from read_temporal_info_to_pcr import getTimedPCRData
from functions_economical_limit_groundwater_extraction import *
from zonalStatistics import zonal_statistics_pcr
from ncRecipes_fixed import getNCDates, createNetCDF

####################
# global variables #
####################

# logger
logger = logging.getLogger(__name__)

# single year spinup: does what it says to speed processing
fixed_yearly_increment = 1

#############
# functions #
#############

def recast_real_as_natural_ratio(r_f):
    '''Returns the ratio r_n which is the ratio of the natural number n, \
so that r_n= 1:n or r_n= n:1 on the basis of the float r_f provided that \
r_f is non zero, otherwise r_n= 0 is retunred.

At the moment, no check on precision is included.

    Input:
    ======
    r_f:                    float

    Output:
    =======
    r_n:                    float, the ratio of 1:n or n:1, or zero,
                            dependent on the sign and value of r_f
    n_s:                    integer, the selected root to which the decimal
                            part of r_f is scaled; n_s * r_f should again
                            return an integer, including zero.

'''

    #-Initialization
    #-check on r_f
    if not isinstance(r_f, float) and not isinstance(r_f, int):
        sys.exit('r_n is not a float or cannot be recast as a float')
    elif isinstance(r_f, int):
        r_f= float(r_f)
    #-initialize r_n and n_s
    # r_n is set to zero and this value is returned if r_f is zero
    n_s= 1
    r_n= 0.0
    # a precision is defined, a value of zero requires an exact match
    # so that r_f * n_s % 1 equals zero
    precision= 1.0e-5
    n_max= int(1 // precision + 1)
    
    #-Start
    # If r_f is not zero, split the absolute number abs(r_f) into its inte-
    # ger root and decimal fraction and keep its sign, test on d and return
    # the new value of r_n
    if r_f != 0:
        #-sign, integer and decimal part of r_f
        if r_f < 0:
            s_f= -1.0
        else:
            s_f=  1.0
        m_r= int(abs(r_f))
        d_r= abs(r_f) % 1
        if d_r <= (precision**0.5 * m_r) and m_r > 0:
            d_r= 0.0
                
        #-recast d_r as n_r
        # Get the absolute inverse of d_r, n_r, and create a list of integers
        # up to the order of n_r + 1, n_list; next, check when the ratio of 
        # an element of n_list, n_s over d_r returns an integer (or a float
        # close to one), then return n_s and n_r to compute the ratio that
        # constitutes the decimal part of r_n= m_r + n_r / n_s
        if d_r > 0:
            # get the list of integers to process
            n_list= list(range(2, n_max+1))
            # start at n_s is 1 and compute n_r,
            # next check whether n_r is an integer and if not, repeat
            # by getting n_s from n_list in ascending order until this
            # condition is met
            n_s= 1
            n_r= d_r * n_s
            while abs(n_r - round(n_r)) > precision and len(n_list) > 0:
                n_s= n_list.pop(0)
                n_r= d_r * n_s
        #-return r_n, the sum of the integer and the decimal part, the 
        # latter being n_r / n_s; the total sum is multiplied by the sign.
        n_r= round(d_r * n_s)
        r_n= s_f* (m_r + n_r / n_s)
    
    #-Return
    #-return r_n and n_s
    return r_n, n_s

def retrieve_ids_from_mask(mask, id_map, id_list, id_names, message_str = ''):

    '''returns the list of selected ids as well as a message string that lists the names'''

    # initialize the selected ids, then check on the mask to see if further
    # processing is necessary
    selected_ids = []

    # this speeds up things:
    if pcr.cellvalue(pcr.mapmaximum(pcr.scalar(mask)), 1)[0] > 0:
        selection_list = zonal_statistics_pcr(mask, \
            id_map, id_list, np.max)
        for selected_id in id_list:
            if selection_list[id_list.index(selected_id)] == 1:
                selected_ids.append(selected_id)
                message_str = str.join('\n', ( \
                    message_str,
                    '\t%3d: %s' % (selected_id, id_names[selected_id]), \
                    ))
    # return selected_ids and message_str
    return selected_ids, message_str

####################
# class definition #
####################

class pcr_economic_limit_groundwater_extraction(pcrm.StaticModel):
    
    #-definition of the static model script object to allocate areas to crops in
    # order to satisfy demand

    ###################
    # binding section #
    ###################
    def __init__(self, modelconfiguration):
        pcrm.StaticModel.__init__(self)

        ##################
        # Initialization #
        ##################
        # Initialization: most information is supplied from the configuration file

        # Parameters
        # dummy variable name and a list of julian days of the first day of the month;
        # this list is used to read in the crop factors
        self.dummyvariablename = 'dummy'
        self.missing_value = -999.9
        self.full_report = modelconfiguration.general['full_report']
        self.testverbose = modelconfiguration.general['testverbose']
        # this is not going to work but it gives results
        self.date_selection_method = 'exact'
        # conversion from tonnes to kg:
        # conversion for producer price: price is in $ per metric tonne, whereas
        # production as computed in kg per cell area
        self.convfactor_tonnes_to_kg = 1000.0
        self.convfactor_ha_to_m2 = None

        # set the pumping efficiency
        self.initial_pumping_efficiency = 0.95
        self.final_pumping_efficiency   = 0.60

        # set the years to process; config info read as list
        self.startyear = int(modelconfiguration.general['startyear'])
        self.max_number_years = int(modelconfiguration.general['max_number_years'])

        # dynamic crop info?
        if 'true' in modelconfiguration.general['dynamic_crop_info_present'].lower():
                self.dynamic_crop_info_present = True
        else:
                self.dynamic_crop_info_present = False
        
        # and a dummy date to read in values
        dummydate = datetime.datetime(self.startyear, 1, 1)

        # setting on reconstruction of existing wells
        self.reconstruct_existing_wells = \
                modelconfiguration.general['reconstruct_existing_wells']

        # Set the clone on the basis of the spatial attributes of the clone map
        self.cloneattributes= spatialAttributes(modelconfiguration.general['clone'])
        setClone(self.cloneattributes)

        # read in the clone map
        self.landmask = pcr.cover(pcr.readmap(modelconfiguration.general['clone']) != 0, \
                                  pcr.boolean(0))

        # read in the countries, regions and cell area
        # countries
        self.countries = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
            modelconfiguration.general['countries'], 'INT32', 'Nominal',\
            self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
            self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
            pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        self.countries = pcr.ifthen(self.landmask, self.countries)
        self.countries = pcr.ifthen(self.countries != 0, self.countries)
        # cell area

        # get the multiplication factor
        dataattributes = spatialAttributes(modelconfiguration.general['cellarea'])
        productResolution, data_scale_division  = recast_real_as_natural_ratio(\
                                                        dataattributes.xResolution)
        productResolution, clone_scale_division = recast_real_as_natural_ratio(\
                                                        self.cloneattributes.xResolution)
        scale_factor = (data_scale_division / clone_scale_division)**2

        # get the cell area
        self.cellarea = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
            modelconfiguration.general['cellarea'], 'FLOAT32', 'Scalar',\
            self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
            self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
            pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        self.cellarea = pcr.ifthen(self.landmask, scale_factor * self.cellarea)

        # regions
        self.regions =  pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
            modelconfiguration.general['regions'], 'INT32', 'Nominal',\
            self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
            self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
            pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        self.regions = pcr.ifthen(self.landmask, self.regions)

        # read in the groundwater information
        # specific yields of the confined and unconfined layers
        # TODO: replace file access by file_handler's read_file_entry method
        #       but this currently is not compatible with python 2 due to the netCDF
        #       class.
        # read in the specific yield
        # unconfined
        if os.path.isfile(modelconfiguration.groundwater['specific_yield_unconflayer']):
           specificyield_unconfined = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
                modelconfiguration.groundwater['specific_yield_unconflayer'], 'FLOAT32', 'Scalar',\
                self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
                self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
                pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        else:
            specificyield_unconfined = modelconfiguration.convert_string_to_input( \
                                            modelconfiguration.groundwater['specific_yield_unconflayer'], float)
        # confined
        if os.path.isfile(modelconfiguration.groundwater['specific_yield_conflayer']):
           specificyield_confined = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
                modelconfiguration.groundwater['specific_yield_conflayer'], 'FLOAT32', 'Scalar',\
                self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
                self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
                pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        else:
            specificyield_confined = modelconfiguration.convert_string_to_input( \
                                            modelconfiguration.groundwater['specific_yield_conflayer'], float)
        # and just make sure they are spatial scalar PCRaster maps
        specificyield_unconfined = pcr.spatial(pcr.scalar(specificyield_unconfined))
        specificyield_confined   = pcr.spatial(pcr.scalar(specificyield_confined))
        # and put them in the list for use
        self.specificyields      = [specificyield_confined, specificyield_unconfined]

        # read in the confined layer depth, the unconfined depth is set to a missing value
        # unconfined
        layerdepth_unconfined = self.missing_value
        # confined
        if os.path.isfile(modelconfiguration.groundwater['confinining_layer_thickness']):
            layerdepth_confined = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
                modelconfiguration.groundwater['confinining_layer_thickness'], 'FLOAT32', 'Scalar',\
                self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
                self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
                pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        else:
            layerdepth_confined = modelconfiguration.convert_string_to_input( \
                                            modelconfiguration.groundwater['confinining_layer_thickness'], float)
        # and just make sure they are spatial scalar PCRaster maps
        layerdepth_confined   = pcr.spatial(pcr.scalar(layerdepth_confined))
        self.layerdepths = [layerdepth_confined, layerdepth_unconfined]

        # odd variables
        self.vars_in_list = {'layerdepth_confined'      : ('layerdepths',     0), \
                             'layerdepth_unconfined'    : ('layerdepths',     1), \
                             'specificyield_confined'   : ('specificyields',  0), \
                             'specificyield_unconfined' : ('specificyields',  1), \
                            }

        # read in the original groundwater depth
        self.groundwaterdepth = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
            modelconfiguration.groundwater['groundwaterdepth'], 'FLOAT32', 'Scalar',\
            self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
            self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
            pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
        self.groundwaterdepth = pcr.ifthen(self.landmask, pcr.max(0, self.groundwaterdepth))
        self.original_groundwaterdepth = self.groundwaterdepth
        self.projected_groundwaterdepth = self.groundwaterdepth

        # read in the output path; creat and set it if necessary
        self.outputpath = modelconfiguration.general['outputpath']
        
        # if it does not exist, create it
        if not os.path.isdir(self.outputpath):
            os.makedirs(self.outputpath)

        # crop information: set the necesary information, then add the completed
        # file names to the dictionaries
        self.crop_systems     = modelconfiguration.convert_string_to_input( \
                                modelconfiguration.crops['crop_systems'], str)
        self.selected_crops   = modelconfiguration.convert_string_to_input( \
                                modelconfiguration.crops['selected_crop_ids'], int)
        self.crop_ids         = modelconfiguration.convert_string_to_input( \
                                modelconfiguration.crops['crop_ids'], int)
        self.crop_names       = modelconfiguration.convert_string_to_input( \
                                modelconfiguration.crops['crop_names'], str)
        self.crop_types       = dict(zip(self.crop_ids, self.crop_names))
        self.irrigation_types = modelconfiguration.convert_string_to_input( \
                                modelconfiguration.crops['irrigation_types'], str)

        # the following variables are initialized als local variables only and converted
        # into class attributes as cell averages
        irrigation_efficiency = dict(zip(self.irrigation_types, \
                                        modelconfiguration.convert_string_to_input( \
                                        modelconfiguration.crops['irrigation_efficiency'], float)))
        irrigation_overpressure_ratio = \
                                dict(zip(self.irrigation_types, \
                                        modelconfiguration.convert_string_to_input( \
                                        modelconfiguration.crops['irrigation_overpressure_ratio'], float)))

        # set cell-averaged characteristics
        # initialize the overall overpressure and the irrigation efficiency
        self.fraction_irrigation_type = {}
        self.head_overpressure = pcr.scalar(0)
        self.irrigation_efficiency = pcr.scalar(0)

        # read in the irrigation fractions and set the efficiency as well as the over pressure
        # overpressure is in Pa and is converted to m
        sprinkler_radius     = modelconfiguration.convert_string_to_input( \
                               modelconfiguration.crops['sprinkler_radius'], float)
        sprinkler_correction = modelconfiguration.convert_string_to_input( \
                               modelconfiguration.crops['sprinkler_correction'], float)
        head_overpressure, velocity = compute_overpressure_sprinklers( \
                                sprinkler_radius, cor_factor = sprinkler_correction)
        head_overpressure = head_overpressure / 1000 / 9.81

        # iterate over the irrigation types
        for irrigation_type in self.irrigation_types:
            
            # get the fraction
            self.fraction_irrigation_type[irrigation_type] = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
                modelconfiguration.crops['fraction_irrigation_type_fileroot'] % irrigation_type, \
                'FLOAT32', 'Scalar',\
                self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
                self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
                pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)        

        # correct the fraction irrigation type to unity where the values are not properly set
        self.fraction_irrigation_type['surface'] = pcr.max(0.0, 1.0 - \
             pcr.cover(self.fraction_irrigation_type['drip'] + \
                       self.fraction_irrigation_type['sprinkler'], 0))

        # iterate over the irrigation types
        for irrigation_type in self.irrigation_types:
            
            # update the overpressure
            self.head_overpressure = self.head_overpressure + \
                    self.fraction_irrigation_type[irrigation_type] * \
                    irrigation_overpressure_ratio[irrigation_type] * head_overpressure
           
            # update the efficiency
            self.irrigation_efficiency = self.irrigation_efficiency + \
                    self.fraction_irrigation_type[irrigation_type] * \
                    irrigation_efficiency[irrigation_type]

        # complement the file names
        crop_water_requirement_fileroot = modelconfiguration.crops['crop_water_requirement_fileroot']
        maximum_crop_area_fileroot = modelconfiguration.crops['maximum_crop_area_fileroot']
        # and initialize the dictionaries
        self.maximum_crop_area_files = {}
        self.crop_water_requirement_files = {}
        # iterate over the crop systems
        for crop_system in self.crop_systems:

            # nest the dictionary for crop systems
            self.maximum_crop_area_files[crop_system] = {}
            self.crop_water_requirement_files[crop_system] = {}
            
            # iterate over the crop types
            for crop_id in self.selected_crops:
                
                # get the name
                crop_name = self.crop_types[crop_id]

                # add the files
                self.maximum_crop_area_files[crop_system][crop_name] = \
                    maximum_crop_area_fileroot % (crop_system, crop_name)
                self.crop_water_requirement_files[crop_system][crop_name] = \
                    crop_water_requirement_fileroot % (crop_system, crop_name)

        # table information to be read from text files: all relevant model
        # configuration items are set here as internal variables so that they
        # can be updated in the initial / dynamic sections of the model
        
        # water productivity
        self.water_productivity_tbl_filename = modelconfiguration.crops['water_productivity_tbl']

        # country information: table with various information organized per coun-
        # try and the entries provided per row
        # country_table_column_start: this is the last column of those columns identifying
        # the different information on country and region IDs as contained by the country-specfic
        # tables:
        # FID   ISO1    ISO2    ISO3    Name    IMAGE26 IMAGE26_name    fraction_agricultural_population
        self.country_data_tbl_filename = modelconfiguration.table_info['country_data_tbl_filename']
        #~ self.agricultural_population_tbl_filename = modelconfiguration.countrysettings['agricultural_population_tbl_filename']

        self.table_info = {}
        for key, value in modelconfiguration.table_info.items():
            if key[:6] == 'table_':
                self.table_info[key[6:]] = int(value)

        #~ self.agricultural_population_tbl_columns = modelconfiguration.countrysettings['agricultural_population_tbl_columns']

        # adjustment info for the sensitivity analysis
        self.adjustment_info = {}
        if 'adjustments' in vars(modelconfiguration).keys():
            
            for key in modelconfiguration.adjustments.keys():
            
                # check on the value
                if modelconfiguration.adjustments[key][:2] == '\-':
                    modelconfiguration.adjustments[key] = modelconfiguration.adjustments[key][1:]

                # process the value
                if os.path.isfile(modelconfiguration.adjustments[key]):
                    
                    # value is a PCRaster map
                    value = pcr.cover(getattr(spatialDataSet(self.dummyvariablename,\
                            modelconfiguration.adjustments[key], 'FLOAT32', 'Scalar',\
                            self.cloneattributes.xLL, self.cloneattributes.xUR, self.cloneattributes.yLL, self.cloneattributes.yUR,\
                            self.cloneattributes.xResolution, self.cloneattributes.yResolution,\
                            pixels= self.cloneattributes.numberCols, lines= self.cloneattributes.numberRows), self.dummyvariablename), 0)
                else:
                
                    value = modelconfiguration.convert_string_to_input( \
                            modelconfiguration.adjustments[key], float)

                self.adjustment_info[key] = value

        # ******************
        # * model products *
        # ******************
        # define the model products, set the settings to create the netCDF files
        # and initialize the netCDF files

        # *********************
        # * netCDF attributes *
        # *********************
        self.ncattributes= {}
        excludedkeys= []
        for key, value in modelconfiguration.netcdfattrs.items():
            if key not in excludedkeys:
                self.ncattributes[key]= value
        if not 'history' in self.ncattributes.keys():
            self.ncattributes['history']= ''
        self.ncattributes['history']+= '\ncreated on %s.' % (datetime.datetime.now()) 

        # resolution, number of rows and columns
        productResolution, scale_division = recast_real_as_natural_ratio(\
                self.cloneattributes.xResolution)
        number_rows= (self.cloneattributes.yUR - self.cloneattributes.yLL + \
            0.5 * productResolution) // productResolution
        number_cols= (self.cloneattributes.xUR - self.cloneattributes.xLL + \
            0.5 * productResolution) // productResolution

        # set latitudes and longitudes
        latitudes =  -np.arange(number_rows) /\
            scale_division + self.cloneattributes.yUR - 0.5 * productResolution
        longitudes =  np.arange(number_cols) /\
            scale_division + self.cloneattributes.xLL + 0.5 * productResolution

        # model products that are written in the approprate places
        # in the dynamic section
        # full_report lists additional variables that can be reported
        self.modelproducts= {}      

        # year groundwater limit is met
        variablename= 'year_groundwater_limit_met'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= '-'
        self.modelproducts[variablename]['method']= 'frequency'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= False
        self.modelproducts[variablename]['report']= True
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # total_extracted_volume
        variablename= 'total_extracted_volume'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'm3'
        self.modelproducts[variablename]['method']= 'frequency'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= True
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # groundwater depth
        variablename= 'groundwaterdepth'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'm'
        self.modelproducts[variablename]['method']= 'frequency'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # number wells
        variablename= 'number_wells'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= '-'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= True
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # well depth
        variablename= 'well_depth'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'm'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= True
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # actual efficiency
        variablename= 'actual_pumping_efficiency'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= '-'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # projected efficiency
        variablename= 'projected_pumping_efficiency'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= '-'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # investment cost per year
        variablename= 'investment_cost'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # energy cost per year
        variablename= 'energy_cost'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$/year'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # interest cost per year
        variablename= 'interest_cost_investment'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$/year'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # total_revenue_loans
        variablename= 'total_revenue_loans'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$/year'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= self.full_report
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # total_revenue_savings
        variablename= 'total_revenue_savings'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$/year'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= False
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # net_present_value_loans
        variablename= 'net_present_value_loans'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$/year'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= True
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # net_present_value_savings
        variablename= 'net_present_value_savings'
        self.modelproducts[variablename]= {}
        self.modelproducts[variablename]['alias']= variablename
        self.modelproducts[variablename]['conversionfactor']= 1.00
        self.modelproducts[variablename]['unit']= 'US$'
        self.modelproducts[variablename]['method']= 'average'
        self.modelproducts[variablename]['data']= variablename
        self.modelproducts[variablename]['dynamic']= True
        self.modelproducts[variablename]['report']= False
        self.modelproducts[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # crop information, this sets fixed data to files
        # in the case of irrigated crops: 
        # - irrigated_crop_area
        # - total_water_requirement_irrigation
        # - total_income_irrigated
        
        # this information is not written by default!
        self.write_crop_info = self.full_report
        self.cropinfo_products = {}
        
        # irrigated crop area
        variablename= 'irrigated_crop_area'
        self.cropinfo_products[variablename]= {}
        self.cropinfo_products[variablename]['alias']= variablename
        self.cropinfo_products[variablename]['conversionfactor']= 1.00
        self.cropinfo_products[variablename]['unit']= 'm2'
        self.cropinfo_products[variablename]['method']= 'average'
        self.cropinfo_products[variablename]['data']= variablename
        self.cropinfo_products[variablename]['dynamic']= True
        self.cropinfo_products[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # total_income_irrigated
        variablename= 'total_income_irrigated'
        self.cropinfo_products[variablename]= {}
        self.cropinfo_products[variablename]['alias']= variablename
        self.cropinfo_products[variablename]['conversionfactor']= 1.00
        self.cropinfo_products[variablename]['unit']= 'US$/year'
        self.cropinfo_products[variablename]['method']= 'average'
        self.cropinfo_products[variablename]['data']= variablename
        self.cropinfo_products[variablename]['dynamic']= True
        self.cropinfo_products[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # total_water_requirement_irrigation
        variablename= 'total_water_requirement_irrigation'
        self.cropinfo_products[variablename]= {}
        self.cropinfo_products[variablename]['alias']= variablename
        self.cropinfo_products[variablename]['conversionfactor']= 1.00
        self.cropinfo_products[variablename]['unit']= 'm3/year'
        self.cropinfo_products[variablename]['method']= 'average'
        self.cropinfo_products[variablename]['data']= variablename
        self.cropinfo_products[variablename]['dynamic']= True
        self.cropinfo_products[variablename]['filename']= os.path.join(\
            self.outputpath, 'netcdf', '%s.nc' % variablename.lower())

        # initialize the netCDF files
        for variablename in self.modelproducts.keys():
            
            if self.modelproducts[variablename]['report']:
                variablelongname= variablename
                createNetCDF(self.modelproducts[variablename]['filename'], \
                    longitudes, latitudes, 'longitude', 'latitude', 'time',\
                    self.modelproducts[variablename]['alias'], self.modelproducts[variablename]['unit'],\
                    self.missing_value, self.ncattributes, varLongName= variablelongname)
                
        if self.write_crop_info:
            for variablename in self.cropinfo_products.keys():

                variablelongname= variablename
                createNetCDF(self.cropinfo_products[variablename]['filename'], \
                    longitudes, latitudes, 'longitude', 'latitude', 'time',\
                    self.cropinfo_products[variablename]['alias'], self.cropinfo_products[variablename]['unit'],\
                    self.missing_value, self.ncattributes, varLongName= variablelongname)

    ######################
    # internal functions #
    ######################
    
    def obtain_crop_info(self, \
            crop_water_requirement_files, \
            maximum_crop_area_files, \
            nc_date_info, date, date_selection_method, \
            water_productivity, \
            crop_systems, selected_crops, crop_types, \
            landmask, cloneattributes, irrigation_efficiency = 0.7):
        
        '''
obtain_crop_info: returns crop information using information on the files and \
dates and the water productivity.

'''

        # obtain the production per crop type as mass per cell
        rainfed_crop_production_info   = {}
        irrigated_crop_production_info = {}

        # set the totals of production and area
        rainfed_crop_area   = pcr.scalar(0)
        irrigated_crop_area = pcr.scalar(0)
        total_crop_area= pcr.scalar(0)
        rainfed_crop_production = pcr.scalar(0)
        irrigated_crop_production = pcr.scalar(0)
        total_crop_production = pcr.scalar(0)

        # total water requirement: only for irrigation
        total_water_requirement_irrigation = pcr.scalar(0)

        # iterate over all crop systems and the 
        # crop types
        for crop_system in crop_systems:

            # iterate over the crop types
            for crop_id in selected_crops:
                
                # get the name
                crop_name = crop_types[crop_id]
        
                # netCDF file name and dates: crop area
                nc_filename = maximum_crop_area_files[crop_system][crop_name]
                nc_date_list = nc_date_info[nc_filename]
                crop_area, returned_date, message_str= \
                    getTimedPCRData(nc_filename,\
                        nc_date_list, date,\
                        dateSelectionMethod= date_selection_method,\
                        dataAttributes= cloneattributes)

                # get the multiplication factor
                dataattributes = spatialAttributes(nc_filename)
                productResolution, data_scale_division  = recast_real_as_natural_ratio(\
                                                                dataattributes.xResolution)
                productResolution, clone_scale_division = recast_real_as_natural_ratio(\
                                                                self.cloneattributes.xResolution)
                scale_factor = (data_scale_division / clone_scale_division)**2

                crop_area = scale_factor * crop_area
              
                # crop areas have been updated on velocity and snellius!
                # ~ # *** temporary ***
                # ~ # correction for crop area
                # ~ crop_area = 100.0 * crop_area
                # ~ # *** temporary ***
              
                # netCDF file name and dates: crop water requirement
                nc_filename = crop_water_requirement_files[crop_system][crop_name]
                nc_date_list = nc_date_info[nc_filename]
                crop_water_requirement, returned_date, message_str= \
                    getTimedPCRData(nc_filename,\
                        nc_date_list, date,\
                        dateSelectionMethod= date_selection_method,\
                        dataAttributes= cloneattributes)

                # get the crop yield
                crop_yield = pcr.ifthen(landmask, \
                    water_productivity[crop_system][crop_id] * \
                    crop_water_requirement)

                # add the data to the dictionary
                if  crop_system == 'rainfed': 
                    rainfed_crop_production_info[crop_id] = crop_yield * crop_area
                    rainfed_crop_area = rainfed_crop_area + crop_area
                    rainfed_crop_production = rainfed_crop_production + \
                        rainfed_crop_production_info[crop_id]

                if  crop_system == 'irrigated':
                    irrigated_crop_production_info[crop_id] = crop_yield * crop_area
                    irrigated_crop_area = irrigated_crop_area + crop_area
                    irrigated_crop_production = irrigated_crop_production + \
                        irrigated_crop_production_info[crop_id]
                    total_water_requirement_irrigation = \
                        total_water_requirement_irrigation + \
                        pcr.ifthen(landmask, \
                            crop_area * crop_water_requirement / irrigation_efficiency)

                    #~ # *** temporary ***
                    #~ # as a test on the data, report the yield in tonnes per ha for
                    #~ # irrigated wheat
                    #~ if crop_id == 1:
                        #~ pcr.report( \
                            #~ 10.0 * irrigated_crop_production_info[crop_id] / crop_area,\
                            #~ 'temp_wheat_yield.map')
                    #~ # *** temporary ***

        # and obtain the total crop area and total production
        total_crop_area = rainfed_crop_area + irrigated_crop_area
        total_crop_production = rainfed_crop_production + irrigated_crop_production

        # return the output
        return rainfed_crop_production_info, rainfed_crop_area, \
                rainfed_crop_production, \
            irrigated_crop_production_info, irrigated_crop_area, \
                irrigated_crop_production, total_water_requirement_irrigation, \
            total_crop_area, total_crop_production

    # end of additional internal functions

    ###################
    # initial section #
    ###################
    def initial(self):

        # echo to screen
        message_str = str.join('\n', \
            (' * Initializing all information that is assumed to be static over time.', \
            'currently this includes the following variables that may become dynamic:', \
            ' - water productivity', \
            ' - producer price per crop', \
            ' - energy costs', \
            ' - labour costs', \
            ' - material costs invariant with depth', \
            ' - material costs per depth', \
            ' - population', \
            ))
        logger.info(message_str)

        # *** truly static information ***
    
        # file information
        # get the dates for all netCDFs
        nc_date_info = {}
        
        # get the file names, and add its dates
        for crop_system in self.crop_systems:

            # iterate over the crop types
            for crop_id in self.selected_crops:
                
                # get the name
                crop_name = self.crop_types[crop_id]
        
                # netCDF file name
                # crop water requirement
                nc_filename = self.crop_water_requirement_files[crop_system][crop_name]
                nc_date_info[nc_filename]= getNCDates(nc_filename)
                # maximum crop area
                nc_filename = self.maximum_crop_area_files[crop_system][crop_name]
                nc_date_info[nc_filename]= getNCDates(nc_filename)
    
        # initializion of internal variables
        # set the amount of money saved to zero as it can be updated during the
        # spinup period, which is equal to the well life time and the present value
        self.total_revenue_loans   = pcr.scalar(0)
        self.total_revenue_savings = pcr.scalar(0)
        self.net_present_value_loans   = pcr.scalar(0)
        self.net_present_value_savings = pcr.scalar(0)

        # set the long-term water requirements
        self.longterm_water_requirement_irrigation = pcr.scalar(0)
        
        # set the total extracted groundwater volume
        self.total_extracted_volume = pcr.scalar(0)

        # initialize the costs at zero
        self.pump_installation_cost = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.well_construction_cost = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.irrigation_installation_cost = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.investment_cost = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.interest_cost_investment = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.energy_cost = pcr.ifthen(self.landmask, pcr.scalar(0))
        
        # initialize current investment level at zero and initialize 
        # the investment year       
        self.current_investment = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.investment_year = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.abstraction_year = pcr.ifthen(self.landmask, pcr.scalar(0))
        
        # set the initial number of pumps to zero and the well depth as well
        self.number_wells = pcr.ifthen(self.landmask, pcr.scalar(0))
        self.well_depth = pcr.ifthen(self.landmask, pcr.scalar(0))
        
        # *** information that may become dynamic ***
        
        # water productivity: dictionary per crop for each crop system
        self.water_productivity = {}
        self.water_productivity['rainfed'] = read_water_productivity( \
                datafile = self.water_productivity_tbl_filename, \
                selectedcrops = self.selected_crops, testverbose = self.testverbose)
        self.water_productivity['irrigated'] = read_water_productivity( \
                datafile = self.water_productivity_tbl_filename, \
                selectedcrops = self.selected_crops, testverbose = self.testverbose)

        # producer price
        self.rainfed_crop_price_info =  {}
        self.irrigated_crop_price_info= {}

        selected_keys = list(self.table_info.keys())

        # iterate over the crops to get the yield and price
        for crop_id in self.selected_crops:
            # get the column index
            data_key    = 'price_%s' % self.crop_types[crop_id]
            data_column = self.table_info[data_key]
            selected_keys.remove(data_key)
            # add the information to the crop information
            self.rainfed_crop_price_info[crop_id] = \
                map_table_info_to_pcr(\
                    datafile      = self.country_data_tbl_filename, \
                    key_map       = self.countries, \
                    key_column    = self.table_info['country_id'], \
                    data_column   = data_column, \
                    pcr_data_type = pcr.scalar, \
                    testverbose   = False)
            self.irrigated_crop_price_info[crop_id] = \
                map_table_info_to_pcr(\
                    datafile      = self.country_data_tbl_filename, \
                    key_map       = self.countries, \
                    key_column    = self.table_info['country_id'], \
                    data_column   = data_column, \
                    pcr_data_type = pcr.scalar, \
                    testverbose   = False)
        # echo to screen
        message_str = 'currently the water productivity is fixed for all years.'
        logger.warning(message_str)

        # mapped info: countries (test) and mapped information that is not crop dependent
        # as stored in the selected keys
        selected_keys.remove('country_name')     
        self.mapped_country_info = {}
 
        # add the country-specific information
        for data_key in selected_keys:
            # get the mapped data
            data_column = self.table_info[data_key]
            self.mapped_country_info[data_key] = map_table_info_to_pcr(\
                    datafile      = self.country_data_tbl_filename, \
                    key_map       = self.countries, \
                    key_column    = self.table_info['country_id'], \
                    data_column   = data_column, \
                    pcr_data_type = pcr.scalar, \
                    testverbose   = False)

        #~ # process the information on the fractional agricultural population   
        #~ for variable_name, data_column in self.agricultural_population_tbl_columns.items():
            #~ # get the mapped data
            #~ self.mapped_country_info[variable_name] =  map_table_info_to_pcr(\
                    #~ datafile = self.agricultural_population_tbl_filename, \
                    #~ key_map = self.countries, \
                    #~ key_column = self.table_country_key, \
                    #~ data_column = data_column, \
                    #~ pcr_data_type = pcr.scalar, \
                    #~ testverbose = False)

        # update the interest rates to fractions 
        for key in ['discount_rate', 'interest_rate_loans', 'interest_rate_savings']:
            self.mapped_country_info[key] = \
                    0.01 * self.mapped_country_info[key]
        
        # this has been added to the table and has been disabled
        ## *** temporary! ***
        ## add savings fraction and interest and discount rate
        #self.mapped_country_info['discount_rate'] = \
        #        pcr.max(0.01, self.mapped_country_info['interest_rate_loans'] - 0.02)
        #self.mapped_country_info['interest_rate_loans'] = self.mapped_country_info['discount_rate'] + 0.02
        #self.mapped_country_info['interest_rate_savings'] = \
        #        self.mapped_country_info['interest_rate_loans']
        #self.mapped_country_info['fraction_saved'] = \
        #        pcr.ifthen(self.landmask, pcr.scalar(1.00))
        ## update the pumping rate to the lower number
        #self.mapped_country_info['maximum_well_capacity'] = \
        #    pcr.spatial(pcr.scalar(770921))  
        ## *** temporary! ***
        
        # get the country names and IDs
        self.country_names = read_dictionary_from_table( \
                    datafile    = self.country_data_tbl_filename, \
                    key_column  = self.table_info['country_id'], \
                    data_column = self.table_info['country_name'], \
                    testverbose = False)
        self.country_ids = list(self.country_names.keys())
        self.country_ids.sort()
        
        # investment status: identifies if institutional loans are available to
        # farmers if True; if False, savings have to be used
        self.investment_status = pcr.boolean(self.mapped_country_info['investment_status'])
        del self.mapped_country_info['investment_status']
    
        logger.warning ('Investment status is uniformly set to %d' % \
            pcr.cellvalue(pcr.mapmaximum(pcr.scalar(self.investment_status)), 1)[0])
        
        # population information
        self.total_food_costs_irrigated = pcr.ifthen(self.landmask, pcr.scalar(0))
        logger.warning('Population information is not read yet!')

        # finally, add the costs of irrigation per m2
        logger.info('Adding costs of irrigation')
        # initialize the total
        self.mapped_country_info['cost_irrigation'] = pcr.scalar(0)
        # iterate over the irrigation types; remove the maps of cost per type
        for irrigation_type in self.irrigation_types:
            
            # get the key
            cost_key = 'cost_irr_%s' % irrigation_type
            
            # add the cost per country, per type
            self.mapped_country_info['cost_irrigation'] = \
                    self.mapped_country_info['cost_irrigation'] + \
                    self.fraction_irrigation_type[irrigation_type] * \
                    self.mapped_country_info[cost_key]
            
            # remove the cost per type
            del self.mapped_country_info[cost_key]

        ############################
        # adjusting variables with #
        # the information provided #
        ############################
        logger.info('Adjusting any settings as specified')
        for key, value in self.adjustment_info.items():

            # get the function and attribute key
            funckey = key[:3]
            attrkey = key[4:]

            if attrkey in vars(self).keys():

                logger.info('modifying variable %s' % attrkey)

                if funckey == 'mul':

                    setattr(self, attrkey, \
                                  value * getattr(self, attrkey))

                elif funckey == 'add':

                    setattr(self, attrkey, \
                                  value + getattr(self, attrkey))
                else:
                    logger.error('%s cannot be processed yet!' % key)

            elif attrkey in self.mapped_country_info.keys():

                logger.info('modifying variable %s of the mapped country data' % attrkey)

                if funckey == 'mul':

                    self.mapped_country_info[attrkey] = value * \
                                                        self.mapped_country_info[attrkey]

                elif funckey == 'add':

                    self.mapped_country_info[attrkey] = value + \
                                                        self.mapped_country_info[attrkey]

                else:
                    logger.error('%s cannot be processed yet!' % key)

            elif attrkey in self.vars_in_list.keys():
                
                logger.info('modifying variable %s of data organized in list' % attrkey)
                
                varname, varpos = self.vars_in_list[attrkey]
                
                mod_value = getattr(self, varname)[varpos]
                
                if funckey == 'mul':
                    
                    mod_value = mod_value * value

                elif funckey == 'add':

                    mod_value = mod_value + value

                else:
                    logger.error('%s cannot be processed yet!' % key)

                # update the value
                if 'layerdepth' in attrkey:
                    
                    # layer depth to be updated
                    self.layerdepths[varpos] = mod_value
                    
                elif 'specificyield' in attrkey:

                    # specific yield to be updated
                    self.specificyields[varpos] = mod_value
                
                else:
                    
                    logger.error('%s cannot be processed yet!' % attrkey)

            elif 'price_' in attrkey[:6]:

                logger.info('modifying variable %s of data on crop price' % attrkey)

                crop_name = attrkey.replace('price_', '')
                if crop_name in self.crop_names:
                    
                    # get the crop id
                    crop_id = self.crop_ids[self.crop_names.index(crop_name)]
                    
                    if funckey == 'mul':

                        self.rainfed_crop_price_info[crop_id]   = value * \
                                                                  self.rainfed_crop_price_info[crop_id]
                        self.irrigated_crop_price_info[crop_id] = value * \
                                                                  self.irrigated_crop_price_info[crop_id]

                    elif funckey == 'add':
                    
                        self.rainfed_crop_price_info[crop_id]   = value + \
                                                                  self.rainfed_crop_price_info[crop_id]
                        self.irrigated_crop_price_info[crop_id] = value + \
                                                                  self.irrigated_crop_price_info[crop_id]
                    else:
                        
                        logger.error('%s cannot be processed yet!' % attrkey)

            else:
                logger.error('%s cannot be processed yet!' % key)

        #~ sys.exit()

        #################
        # QUASI-DYNAMIC #
        #################

        # *** time management ***
        # create a list of years, that starts prior to the start date by the
        # longest well life time; also set a map of the actual start year of the
        # spinup per cell at which the accumulation of revenue will start.
        
        # initialize the spinup       
        startyear_spinup = self.startyear - self.mapped_country_info['well_lifetime']
        
        # set the first and final year and create a list of years
        firstyear = int(pcr.cellvalue(pcr.mapminimum(startyear_spinup), 1)[0])
        lastyear  = self.startyear + self.max_number_years
        self.years = list(range(firstyear, self.startyear, fixed_yearly_increment)) + \
                     list(range(self.startyear, lastyear + 1, fixed_yearly_increment)) 

        # echo to screen
        message_str = str.join(' ', ( \
            '\n  * The simulation will run over %d years,' % len(self.years), \
            'starting in %d, spinning up to %d and ending in %d' % \
                (firstyear, self.startyear, lastyear)))
        logger.info(message_str)

        # get the weights for the spinup years and actual years
        update_weight_spinup  = pcr.spatial(pcr.scalar( \
            1.0 / (self.startyear - startyear_spinup)))
        update_weight_general = pcr.min(fixed_yearly_increment / self.mapped_country_info['well_lifetime'])
        update_weight_groundwater = pcr.min(fixed_yearly_increment, self.mapped_country_info['well_lifetime'])
        
        # and initialize the general spinup status; once False, all values are
        # updated per year and set a mask that keeps track of the spinup conditions
        spinup_mask = pcr.boolean(0)
        general_spinup_flag = False

        # set the year the groundwater limit is reached; this is set outside the
        # period of interest initially
        self.year_groundwater_limit_met = pcr.ifthen(self.landmask, \
                pcr.scalar(self.years[-1] + 1))
                
        # iterate over the years
        for year in self.years:
            
            # set the date
            date = datetime.datetime(year, 1, 1)
            poscnt = max(0, self.years.index(year) - \
                self.years.index(self.startyear))

            # echo to screen
            message_str = ' * processing %d' % year
            logger.info(message_str)

            ###########################
            # *** time management *** #
            ###########################
            # boolean checks and updates on masks
            # spinup: check the general status and the condition per country
            general_spinup_flag = year < self.startyear
            # echo to screen
            if general_spinup_flag:
                logger.info('   as spinup')
            else:
                logger.info('   as actual year')
            # spinup occurs if the year is between the start year for the particular
            # country and the year less than the start year
            update_spinup_mask = pcr.pcrnot(spinup_mask) & (year >= startyear_spinup) & \
                (year < self.startyear)

            # retrieve selected IDs and the corresponding message str
            message_str = ' * revenue is initialized for the following countries:'            
            selected_ids, message_str = retrieve_ids_from_mask(update_spinup_mask, \
                self.countries, self.country_ids, self.country_names, message_str)
            # echo to screen
            if len(selected_ids) > 0:
                logger.info(message_str)
            # set the spinup mask
            spinup_mask = spinup_mask | update_spinup_mask
            
            ##########################
            # yield, production      #
            # and water requirements #
            ##########################
            # get the information on the production
            # here all information regarding crop production is being processed;
            # this concerns the production, the crop area for both crop systems
            # (rainfed, irrigated) and as total as well as the irrigation water
            # requirement.
          
            # retrieve all crop production info for the current year
            # or the first year only, dependent on the settings
            if year == self.years[0] or self.dynamic_crop_info_present:

                # set the crop information within the dynamic section
                rainfed_crop_production_info, \
                rainfed_crop_area, \
                rainfed_crop_production, \
                irrigated_crop_production_info, \
                irrigated_crop_area, \
                irrigated_crop_production, \
                total_water_requirement_irrigation, \
                total_crop_area, \
                total_crop_production = \
                    self.obtain_crop_info( \
                        self.crop_water_requirement_files, \
                        self.maximum_crop_area_files, \
                        nc_date_info, date, self.date_selection_method, \
                        self.water_productivity, \
                        self.crop_systems, \
                        self.selected_crops, self.crop_types, \
                        self.landmask, self.cloneattributes, \
                        self.irrigation_efficiency, \
                        )

                # all crop information obtained              
                # echo to screen
                message_str = ' - all information on crop production updated for %d' % year
                if not self.dynamic_crop_info_present:
                    message_str = str.join(', ', (message_str, 'and is kept fixed throughout as dynamic crop information is disabled!'))
                logger.info(message_str)

            # update on investments
            # update of well depth if the maximum lifetime of the well is reached
            # and the current investment per year is not exceed the revenue over
            # the past period
            evaluate_welldepth_mask = (pcr.scalar(year) >= self.startyear) & \
                ((pcr.scalar(year - self.startyear) % self.mapped_country_info['well_lifetime']) == 0) & \
                (total_water_requirement_irrigation > 0)

            # retrieve selected IDs and the corresponding message str
            message_str = ' * well depth will be updated for the following countries:'
            selected_ids, message_str = retrieve_ids_from_mask(evaluate_welldepth_mask, \
                self.countries, self.country_ids, self.country_names, message_str)
            # echo to screen
            if len(selected_ids) > 0:
                logger.info(message_str)

            ######################
            # income and revenue #
            ######################

            # get the income
            # here the income is computed for the production of rainfed and
            # irrigated crops. the costs for subsistence can subsequenly be subtract-
            # ed from this, depending on the population that for its livelihood
            # relies on agriculture, to obtain the revenue that can be used to
            # determine the economic limit.

            # echo to screen
            message_str = str.join('\n', (' - updating information on income', \
                '\twarning: currently the producer price per crop is fixed for all years.'))
            logger.info(message_str)

            # compute the total income
            total_income_rainfed = \
                (sum_list(multiply_dicts(rainfed_crop_production_info, \
                    self.rainfed_crop_price_info).values())) / \
                    self.convfactor_tonnes_to_kg
            total_income_irrigated =  \
                (sum_list(multiply_dicts(irrigated_crop_production_info, \
                    self.irrigated_crop_price_info).values())) / \
                    self.convfactor_tonnes_to_kg
            total_income = total_income_rainfed + total_income_irrigated

            # compute the average crop price
            self.producer_price_average = total_income_irrigated / \
                                     sum_list(irrigated_crop_production_info.values()) * \
                                     self.convfactor_tonnes_to_kg

            if self.write_crop_info and not general_spinup_flag:
                for variablename, pcr_field in { \
                        'irrigated_crop_area': irrigated_crop_area, \
                        'total_water_requirement_irrigation': total_water_requirement_irrigation, \
                        'total_income_irrigated': total_income_irrigated, \
                        }.items():
                    # write the output to netCDF
                    writeField(self.cropinfo_products[variablename]['filename'], \
                            pcr.pcr2numpy(pcr.scalar(pcr_field), self.missing_value), \
                                self.cropinfo_products[variablename]['alias'],\
                                date, poscnt, 'time', self.missing_value)
          
            # revenue computed
            logger.info(' - all information on income updated')

            #############################################
            # processing: spinup or general condtition? #
            #############################################
            if general_spinup_flag:

                ##################
                # update revenue #
                # and water      #
                # requirement in #
                # the case of    #
                # spinup         #
                ##################

                # set the spinup year
                spinup_year = year - self.years[0]
                # echo update
                logger.info('   spinup: updating information for %d: %d year in spinup' % \
                    (year, spinup_year))
                # update the longterm revenue and longterm water requirement for those
                # countries where an update is required, i.e., the spinup conditions
                # are True in general and for the country in particular
                self.total_revenue_loans = self.total_revenue_loans + pcr.ifthenelse( \
                        self.investment_status & (spinup_year <= self.mapped_country_info['well_lifetime']), \
                            pcr.max(0, total_income_irrigated - self.total_food_costs_irrigated) / \
                                (1.0 + self.mapped_country_info['discount_rate']) ** spinup_year, \
                            pcr.scalar(0))
                self.total_revenue_savings = self.total_revenue_savings + pcr.ifthenelse( \
                        pcr.pcrnot(self.investment_status) & (spinup_year <= self.mapped_country_info['well_lifetime']), \
                            self.mapped_country_info['fraction_saved'] * \
                                pcr.max(0, total_income_irrigated - self.total_food_costs_irrigated) * \
                                (1.0 + self.mapped_country_info['interest_rate_savings']) ** \
                                (self.mapped_country_info['well_lifetime'] - spinup_year), \
                            pcr.scalar(0))
                #self.total_revenue_savings = total_income_irrigated
                self.longterm_water_requirement_irrigation =  \
                        self.longterm_water_requirement_irrigation + \
                        pcr.ifthenelse(spinup_mask, update_weight_spinup, 0) * \
                            total_water_requirement_irrigation

            else:

                #####################
                # split the update  #
                # mask into areas   #
                # where the invest- #
                # ments are made    #
                # from savings and  #
                # where from loans  #
                #####################

                # echo update
                logger.info('   actual year: updating information for %d' % year)
                
                # get the areas where the updates are to be made from savings
                # and from loans
                # presently, this is updated here on the basis of the initial,
                # fixed map: investment status is True if instutional loans
                # are accessible to farmers.
                # this can be made dynamically here
                investment_status = self.investment_status

                ################################
                # groundwater:                 #
                # actual water depth and       #
                # the projected depth from the #
                # current groundwater depth    #
                # and the corresponding        #
                # pumping efficiency with the  #
                # number of pumps required     #
                ################################

                ################################################################
                ## old computation without confining layer
                ##    
                ## project the groundwater depth in case an update is required
                ##self.projected_groundwaterdepth = pcr.ifthenelse( \
                ##    evaluate_welldepth_mask, \
                ##        self.groundwaterdepth + self.mapped_country_info['well_lifetime'] * \
                ##            self.longterm_water_requirement_irrigation / \
                ##            (self.cellarea * self.porosity), \
                ##        self.projected_groundwaterdepth)
                ################################################################
                
                # project the groundwater depth in case an update is required
                self.projected_groundwaterdepth = pcr.ifthenelse( \
                    evaluate_welldepth_mask, \
                        compute_groundwaterdepth_from_extraction( \
                            self.mapped_country_info['well_lifetime'] * \
                                self.longterm_water_requirement_irrigation, \
                            self.groundwaterdepth, self.cellarea, \
                            self.specificyields, self.layerdepths, mv_id = self.missing_value, \
                            testverbose = False), \
                        self.projected_groundwaterdepth)

                # and in the case of a confining layer, then set the depth at
                # least to the thickness of the confining layer
                self.projected_groundwaterdepth = pcr.ifthenelse( \
                        evaluate_welldepth_mask & (self.layerdepths[0] > 0), \
                            pcr.max(self.layerdepths[0], self.projected_groundwaterdepth), \
                            self.projected_groundwaterdepth)

                # pumping efficency
                self.projected_pumping_efficiency = compute_pumping_efficiency( \
                    self.projected_groundwaterdepth + self.head_overpressure, \
                    self.initial_pumping_efficiency, self.final_pumping_efficiency)
                self.actual_pumping_efficiency = compute_pumping_efficiency( \
                    self.groundwaterdepth + self.head_overpressure, \
                    self.initial_pumping_efficiency, self.final_pumping_efficiency)

                # number of new pumps
                number_new_wells = pcr.ifthenelse( \
                    evaluate_welldepth_mask, \
                    pcr.max(0, self.longterm_water_requirement_irrigation / \
                        (self.projected_pumping_efficiency * \
                        self.mapped_country_info['maximum_well_capacity']) -\
                        self.number_wells), pcr.scalar(0))

                # increase in well depth
                new_well_depth = self.projected_groundwaterdepth + self.mapped_country_info['well_water_depth']

                #####################
                # compute all costs #
                # for the well      #
                #####################
              
                # cost to install the pump at all wells!
                # first, get the pump capacity in MJ per year and convert it to kW
                pump_capacity = compute_total_energy_well_extraction( \
                        self.mapped_country_info['maximum_well_capacity'], \
                        self.projected_groundwaterdepth + self.head_overpressure, \
                        self.projected_pumping_efficiency)
                pump_capacity = pump_capacity  * (1000.0 / (365*86400))
                # get the cost per pump
                self.pump_installation_cost = self.mapped_country_info['pump_cost_start'] + \
                    self.mapped_country_info['pump_cost_kw'] * pump_capacity
                # and scale it to the total number of wells if the update is true
                self.pump_installation_cost = pcr.ifthenelse( \
                        evaluate_welldepth_mask, \
                            self.pump_installation_cost * (self.number_wells + number_new_wells), \
                            pcr.scalar(0))

                # cost to install irrigation for all wells!
                self.irrigation_installation_cost = pcr.ifthenelse( \
                        evaluate_welldepth_mask, \
                            irrigated_crop_area * \
                                self.mapped_country_info['cost_irrigation'], \
                            pcr.scalar(0))               
                
                # energy cost for operation of well per year
                self.energy_cost = compute_specific_energy_cost_withdrawal( \
                        total_water_requirement_irrigation, \
                        self.groundwaterdepth + self.head_overpressure, \
                        self.mapped_country_info['electricity_price'], \
                        self.actual_pumping_efficiency)

                # cost to create the well
                # from existing wells and new wells
                
                # new wells
                new_wells_construction_costs = pcr.ifthenelse( \
                    evaluate_welldepth_mask, \
                        compute_well_costs_with_depth( \
                            self.mapped_country_info['material_cost'], \
                            self.mapped_country_info['labour_cost'], \
                            new_well_depth, number_new_wells), \
                        pcr.scalar(0))
                            
                # existing wells
                if self.reconstruct_existing_wells:
                    existing_wells_reconstruction_costs = pcr.ifthenelse( \
                        evaluate_welldepth_mask, \
                            compute_well_costs_with_depth( \
                                self.mapped_country_info['material_cost'], \
                                self.mapped_country_info['labour_cost'], \
                                new_well_depth, self.number_wells), \
                            pcr.scalar(0))
                else:
                    existing_wells_reconstruction_costs = pcr.ifthenelse( \
                        evaluate_welldepth_mask, \
                            compute_well_costs_with_depth( \
                                self.mapped_country_info['material_cost'], \
                                self.mapped_country_info['labour_cost'], \
                                pcr.max(0.0, new_well_depth - self.well_depth),\
                                self.number_wells), \
                            pcr.scalar(0))

                # total costs of well construction
                self.well_construction_cost = pcr.ifthenelse( \
                    evaluate_welldepth_mask, \
                        existing_wells_reconstruction_costs + \
                            new_wells_construction_costs, \
                        self.well_construction_cost)
                       
                ##########################
                # compute and update     #
                # the investment         #
                # and the annual cost    #
                ##########################
                              
                # current investment set
                self.current_investment = pcr.ifthenelse(evaluate_welldepth_mask, \
                    self.well_construction_cost + self.pump_installation_cost +\
                    self.irrigation_installation_cost, \
                    self.current_investment)
                
                # and the costs per year over the current investment period computed
                self.investment_cost = self.current_investment / \
                    self.mapped_country_info['well_lifetime']

                ############################
                # decide on and            #
                # set the year the limit   #
                # is reached on the basis  #
                # of the revenue and the   #
                # current investment       #
                ############################
                # update the year the well depth on the
                # following conditions:
                # 1) the groundwater limit is not yet reached;
                # 2) the revenue is greater than the investment for the two
                #    selected economic systems (loans, savings);
                # 3) the evaluation of well depth is True.
                # general
                update_welldepth_mask = evaluate_welldepth_mask & \
                    (self.year_groundwater_limit_met >= year)
                # dependent on investment type
                update_welldepth_mask = update_welldepth_mask & \
                    pcr.ifthenelse(investment_status, \
                        self.total_revenue_loans > self.current_investment, \
                        self.total_revenue_savings > self.current_investment)

                ##############################
                # update the number of wells #
                # and their depth            #
                ##############################
                # update the well depth and number of wells
                self.well_depth = pcr.ifthenelse(update_welldepth_mask, \
                        pcr.max(self.well_depth, new_well_depth), self.well_depth)
                self.number_wells = self.number_wells + pcr.ifthenelse(update_welldepth_mask, \
                        number_new_wells, 0)

                ###############################
                # set the year the economical #
                # limit is reached            #
                ###############################
                # this is set as the minimum of the year that is
                # already stored and updated if the well depth is not updated
                # but it is evaluated
                self.year_groundwater_limit_met = pcr.ifthenelse( \
                    evaluate_welldepth_mask & pcr.pcrnot(update_welldepth_mask), \
                        pcr.min(year, self.year_groundwater_limit_met), \
                        self.year_groundwater_limit_met)

                #############################
                # next, update the settings #
                # to continue with the next #
                # round of investments      #
                #############################

                ##############################
                # update the investment year #
                ##############################                
                # reset the investment year
                self.investment_year = pcr.ifthenelse( \
                    self.investment_year == self.mapped_country_info['well_lifetime'], \
                    pcr.scalar(0), self.investment_year)

                # update the investment year
                self.investment_year = self.investment_year + 1
                
                # update the abstraction year
                self.abstraction_year = self.abstraction_year + 1

                ##################################
                # cost of interest on investment #
                ##################################
                self.interest_cost_investment = (1.0 - \
                    (self.investment_year - 1) / self.mapped_country_info['well_lifetime']) * \
                    self.mapped_country_info['interest_rate_loans'] * \
                    self.current_investment

                ################################
                # reset and update the         #
                # total revenue and compute    #
                # the long-term water          #
                # irrigation water requirement #
                ################################
                
                # reset the running mean if the well life time is up
                self.total_revenue_loans = pcr.ifthenelse(evaluate_welldepth_mask, 0, \
                        self.total_revenue_loans)
                self.total_revenue_savings = pcr.ifthenelse(evaluate_welldepth_mask, 0, \
                        self.total_revenue_savings)
                self.longterm_water_requirement_irrigation = \
                        pcr.ifthenelse(evaluate_welldepth_mask, 0, \
                        self.longterm_water_requirement_irrigation)

                # reset the long-term irrigation water requirement:
                # are we still pumping?
                groundwater_pumping_mask = self.year_groundwater_limit_met > year

                # long term irrigation water requirement
                self.longterm_water_requirement_irrigation =  \
                        self.longterm_water_requirement_irrigation + \
                        pcr.ifthenelse(groundwater_pumping_mask, \
                            update_weight_general * total_water_requirement_irrigation, 0)

                # and the revenue to meet the investment via loans or savings
                self.total_revenue_loans = self.total_revenue_loans + pcr.ifthenelse( \
                    self.investment_status & groundwater_pumping_mask, \
                        pcr.max(0, total_income_irrigated - self.total_food_costs_irrigated - \
                            self.interest_cost_investment - self.energy_cost) * \
                            (1.0 + self.mapped_country_info['discount_rate']) ** -self.investment_year, \
                        pcr.scalar(0))
                self.total_revenue_savings = self.total_revenue_savings + pcr.ifthenelse( \
                        pcr.pcrnot(self.investment_status) & groundwater_pumping_mask, \
                            self.mapped_country_info['fraction_saved'] * \
                                 pcr.max(0,total_income_irrigated - self.total_food_costs_irrigated - \
                                 self.energy_cost) * \
                                (1.0 + self.mapped_country_info['interest_rate_savings']) ** \
                                (self.mapped_country_info['well_lifetime'] - self.investment_year), \
                            pcr.scalar(0))
                
                # and make Marc happy with the net present value
                self.net_present_value_loans = self.net_present_value_loans + \
                    pcr.ifthenelse(self.investment_status & groundwater_pumping_mask, \
                        pcr.max(0, total_income_irrigated - self.total_food_costs_irrigated - \
                                self.interest_cost_investment - self.energy_cost - self.investment_cost) * \
                        (1 + self.mapped_country_info['discount_rate'])**(-self.abstraction_year), pcr.scalar(0))
                               
                self.net_present_value_savings = self.net_present_value_savings + \
                    pcr.ifthenelse(pcr.pcrnot(self.investment_status) & groundwater_pumping_mask, \
                        pcr.max(0, total_income_irrigated - self.total_food_costs_irrigated - \
                                self.energy_cost - self.investment_cost) * \
                        (1 + self.mapped_country_info['discount_rate'])**(-self.abstraction_year), pcr.scalar(0))

                # finally, update the groundwater depth if still pumping
                groundwaterdepth_old = self.groundwaterdepth

                # adjust the groundwater depth in case an update is required
                self.groundwaterdepth = pcr.ifthenelse( \
                    groundwater_pumping_mask, \
                        compute_groundwaterdepth_from_extraction( \
                            fixed_yearly_increment * \
                            total_water_requirement_irrigation, \
                            self.groundwaterdepth, self.cellarea, \
                            self.specificyields, self.layerdepths, mv_id = self.missing_value, \
                            testverbose = False), \
                        self.groundwaterdepth)

                # NOTE: this is disabled as the groundwater depth is limited by
                # the extraction and the well depth explicitly considers the
                # thickness of the confining layer
                ## in the case of a confining layer, then set the depth at
                ## least to the thickness of the confining layer
                #self.groundwaterdepth = pcr.ifthenelse( \
                #        groundwater_pumping_mask & (self.layerdepths[0] > 0), \
                #            pcr.max(self.layerdepths[0], self.groundwaterdepth), \
                #            self.groundwaterdepth)

                # and check against the well depth if this is not zero
                self.groundwaterdepth = pcr.ifthenelse(self.well_depth > 0, \
                    pcr.min(self.well_depth, self.groundwaterdepth), self.groundwaterdepth)

                # get the total extracted groundwater volume
                self.total_extracted_volume = self.total_extracted_volume  + \
                    pcr.ifthenelse(groundwater_pumping_mask & \
                            (total_water_requirement_irrigation > 0), \
                            fixed_yearly_increment * \
                            total_water_requirement_irrigation, pcr.scalar(0))

                # and a patch to set the year to the start year if the well depth
                # is still zero and not updated
                if year == self.startyear:
                    logger.info('setting year for areas where no groundwater can be extracted')
                    self.year_groundwater_limit_met = pcr.ifthenelse(\
                                                        self.well_depth == 0, year, \
                                                        self.year_groundwater_limit_met)

                ##################
                # write output   #
                # to netCDF file #
                ##################

                for variablename in self.modelproducts.keys():

                    if self.modelproducts[variablename]['dynamic'] and \
                            self.modelproducts[variablename]['report']:
                        
                        # echo to screen
                        logger.debug(' - writing info on %s to netCDF file for %s' % \
                                (variablename, date))

                        # write to netCDF
                        writeField(self.modelproducts[variablename]['filename'], \
                                pcr.pcr2numpy(pcr.scalar(getattr(self, self.modelproducts[variablename]['data'])), \
                                    self.missing_value), self.modelproducts[variablename]['alias'],\
                                    date, poscnt, 'time', self.missing_value)

                ####################
                # additional tests #
                ####################
                # test on the costs and the projected groundwater depth
                if year == self.startyear and self.testverbose:
                    
                    # report all variables
                    for variablename, value in vars(self).items():
                        if isinstance(value, pcr.Field):
                            pcr.report(value, os.path.join(\
                                       self.outputpath, 'temp', '%s.map' % variablename.lower()))

                    # report the variables for the senstivity analysis
                    for variablename, value in { \
                                'producer_price'             : self.producer_price_average, \
                                'discount_rate'              : self.mapped_country_info['discount_rate'], \
                                'interest_rate_loans'        : self.mapped_country_info['interest_rate_loans'], \
                                'labour_cost'                : self.mapped_country_info['labour_cost'], \
                                'material_cost'              : self.mapped_country_info['material_cost'], \
                                'cost_irrigation'            : self.mapped_country_info['cost_irrigation'], \
                                'electricity_price'          : self.mapped_country_info['electricity_price'], \
                                'layerdepth_confined'        : self.layerdepths[0], \
                                'specificyield_confined'     : self.specificyields[0], \
                                'specificyield_unconfined'   : self.specificyields[1], \
                                'pump_installation_cost'     : self.pump_installation_cost, \
                                'initial_pumping_efficiency' : pcr.scalar(self.initial_pumping_efficiency), \
                                'final_pumping_efficiency'   : pcr.scalar(self.final_pumping_efficiency), \
                                }.items():

                        if isinstance(value, pcr.Field):
                            pcr.report(value, os.path.join(\
                                       self.outputpath, 'temp', '%s.map' % variablename.lower()))

                if self.testverbose:
                    pcr.report(self.year_groundwater_limit_met, os.path.join(\
                               self.outputpath, 'temp', '%s.map' % 'year_groundwater_limit'))

            #~ # quit on first year
            #~ sys.exit('quit on first year for development')

        #############################
        # dynamic section completed #
        #############################

        # any static output, set poscnt to zero
        poscnt = 0
        for variablename in self.modelproducts.keys():

            if not self.modelproducts[variablename]['dynamic'] and \
                        self.modelproducts[variablename]['report']:
                
                # echo to screen
                logger.debug(' - writing info on %s to netCDF file for %s' % \
                        (variablename, date))

                # write to netCDF
                writeField(self.modelproducts[variablename]['filename'], \
                        pcr.pcr2numpy(pcr.scalar(getattr(self, self.modelproducts[variablename]['data'])), \
                            self.missing_value), self.modelproducts[variablename]['alias'],\
                            date, poscnt, 'time', self.missing_value)

        # end of iteration over the years
        message_str = ' * all years processed'
        logger.info(message_str)

#############
# call main #
#############

def main():
    pass

if __name__ == "__main__":
    
    main()

    sys.exit('model ended')
