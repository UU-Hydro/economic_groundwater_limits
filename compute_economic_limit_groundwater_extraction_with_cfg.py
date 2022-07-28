#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
# Modules #
###########
# Modules os and sys are imported by default
#-general modules and packages
import os
import sys
import optparse
import logging

import pcraster as pcr
import pcraster.framework as pcrm

# model specific packages
from model_configuration import configuration_parser

from pcr_economic_limit_groundwater_extraction import pcr_economic_limit_groundwater_extraction

####################
# global variables #
####################

# logger
logger = logging.getLogger(__name__)

# type set to identify None (compatible with pytyon 2.x)
NoneType = type(None)

########
# main #
########

def main():

    #-test specification of configuration file
    usage = 'usage: %prog CFGFILE'
    parser = optparse.OptionParser(usage = usage)
    (options, arguments)= parser.parse_args()
    if len(arguments) < 1:
        parser.error('incorrect number of arguments specified')
    else:   
        cfgfilename = arguments[0]
        subst_args  = arguments[1:]

    cfgfilename = os.path.abspath(cfgfilename)

    # set the configuration object
    # object to handle configuration/ini file
    sections = ['general', 'netcdfattrs', 'groundwater', 'crops', 'table_info']
    groups= []

    model_configuration = configuration_parser(cfgfilename = cfgfilename, \
               sections =  sections, groups = groups, subst_args = subst_args)

    # change to the scratch path
    startdir = os.getcwd()
    os.chdir(model_configuration.temppath)

    # call and run the model to allocate production to the available resources in an area
    run_model= pcrm.StaticFramework(pcr_economic_limit_groundwater_extraction(model_configuration))
    run_model.run()

    # model finished

    # reset directory
    os.chdir(startdir)


if __name__ == '__main__':
    main()
    sys.exit('\nmodel run has finished normally')

