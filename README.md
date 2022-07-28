# github repository "economic_groundwater_limits"
#
# Global Economic Limits of Groundwater used for Irrigation
#
# Copyright (c) L.P.H. (Rens) van Beek / Marc F.P. Bierkens 2018-2022
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
# 
# with questions, please contact r.vanbeek@uu.nl
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
# This repository contains the scripts required to simulate
# whether the economical limit to groundwater exploitation is
# reached or not, conditional to the fact that the total costs for
# exploitation are larger than the revenue.
#
# Included scripts are:
# compute_economic_limit_groundwater_extraction_with_cfg.py
#       script wrapping around the main script handling the input
#       and the model settings by means of the configuration file;
# economic_limit_groundwater_extraction_sensitivity_analysis.cfg
#       an example of the configuration file with wild cards included
#       that can be used as the basis for the senstivity analysis;
#       wildcards are identified by "$1 ... $N" with the number corres-
#       ponding with the position on the command line after the config-
#       uration file;
# pcr_economic_limit_groundwater_extraction.py
#       the main script, parsing the input, spinning up, and executing
#       the simulation for the number of years specified.
#
# Included are the following ancillary scripts:
# functions_economical_limit_groundwater_extraction.py:
#       module with additional functions supporting the analysis;
# basic_functions.py:
#       script providing some additional functions supporting the anal-
#       ysis;
# model_configuration.py:
#       script that parses the information from the configuration file;
# spatialDataSet2PCR.py:
#       script handling spatial input from different file formats, in-
#       cluding netCDF files.
#
# to run the script python 3 is required and several packages, including
# gdal and pcraster. The necessary packages can be installed via the
# yaml file pcraster38.yml included.
#
# To run the model code for the CONUS and surrounding 
# basins, use the following command line with the settings for the most
# likely estimates for the different variables (over multiple lines "\"):

python /myhome/economic_limit_irrigation/compute_economic_limit_groundwater_extraction_with_cfg.py \
        /myhome/economic_limit_irrigation/economic_limit_groundwater_extraction_sensitivity_analysis.cfg \
        "/myhome/economic_limit_irrigation/data/masks/conus_mask.map" \
        "/scratch/fosm_economic_limit/conus_00" "900" \
        "/myhome/economic_limit_irrigation/data/groundwater/confining_layer_thickness_updated.map" \
        "1.0" "1.0" "0.001" "1.0" "0.0" "1.0" "1.0" "1.0" "0.0" "True"

# set up the sensitivity analysis as a test with the following wild cards
# python compute_economic_limit_groundwater_extraction_with_cfg.py $1 $2 ... $N
# $0  : configuration file
# $1  : map of the relevant land mask
# $2  : output directory
# $3  : length of simulation period [years]
# $4  : thickness of the unconfined layer [m]
# $7  : storativity of the confined layer [m/m]
# The following are specific settings for the sensitivity analysis
# $5  : changing the thickness of the unconfined layer;
#       to be decided! either spatially variable or constant.
# $6  : muliplicative factor for the porosity, based on the standard deviation of the global
#     : geological information; symmetric around the mean, spatially variable;
#       it is additive to the log-transformed mean.
# $8  : multiplication factor for the producer prices: applies to the following crops
#       PriceWheat; PriceMaize; PriceRice; PriceBarley; PriceRye; PriceMillet;
#       PriceSorghum; PriceSoybeans; PriceSunflower; PricePotatoes; PriceCassava;
#       PriceSugarCane; PriceSugarbeet; PriceOilpalm; PriceRapeseed; PriceGoundnuts;
#       PricePulses; PriceCitrus; PriceDatepalm; PriceGrapes; PriceCotton; PriceCocoa;
#       PriceCoffee; PriceOthersPerennial; PriceFodderGrasses
#       Note: not all crops may be included; currently others_perennial,
#             managed_grassland_or_pasture, and others_annual are excluded.
#       Uniform multiplication factor, symmetric around the mean and is applied
#       equally to all crops.
# $9  : addition to the discount rate and interest rate:
#       DiscountRate[%]; InterestRateInvestment[%]
#       addition is uniform but asymmetric and applied to the mean.
# $10 : pre-multiplier for the electricity costs: CostElectricity[$/kwh]
#       Uniform multiplication factor, symmetric around the mean.
# $11 : pre-multiplier for the labour costs: CostLabour[$/m]
#       Uniform multiplication factor, symmetric around the mean.
# $12 : pre-multiplier for the material costs:
#       CostMaterials[$/m]; CostPumpStart[$/piece]; CostPumpKW[$/kW]; CostIrrigationSurface [$/m2]; CostIrrigationSprinkler [$/m2]; CostIrrigationDrip [$/m2]
#       Uniform multiplication factor, symmetric around the mean.
# $13 : additional change on the pumping efficiency: PumpingEfficiency[-]
#       addition is uniform and symmetric around the mean.
#
# $14 : boolean variable, if True, additional output files are generated.
#
# end of README.md
