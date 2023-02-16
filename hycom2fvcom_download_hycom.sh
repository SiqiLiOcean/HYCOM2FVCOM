#!/bin/bash

#==========================================================================
# Download HYCOM dataset 
#
# https://www.hycom.org/dataserver
#
# First use explorer to find which dataset covers your reseaching time:
#   https://tds.hycom.org/thredds/catalog.html
#
#
# resolution: 0.08 x 0.08 
# time interval: 3 hours
# format: NetCDF
#
# Require library: NCO (https://nco.sourceforge.net/)
#
# Suppoted GOFS 3.1 dataset:
# 
# | expt  | starting time    | ending time      |
# ----------------------------------------------------
# | v53.X | 1994-01-01 12:00 | 2015-12-31 09:00 |
# | v56.3 | 2014-07-01 12:00 | 2016-09-30 09:00 |
# | v57.2 | 2016-05-01 12:00 | 2017-02-01 09:00 |
# | v92.8 | 2017-02-01 12:00 | 2017-06-01 09:00 |
# | v57.7 | 2017-06-01 12:00 | 2017-10-01 09:00 |
# | v92.9 | 2017-10-01 12:00 | 2018-01-01 09:00 |
# | v93.0 | 2018-01-01 12:00 | 2020-02-19 09:00 |
# | u93.0 | 2018-09-19 12:00 | 2018-12-09 09:00 |
# | y93.0 | 2018-12-04 12:00 | present          |
#
#
# Siqi Li, SMAST
# 2022-11-30
#
# Updates
# 2022-12-16  Siqi Li  Add the settings of longitude/latitude limits
# 2023-01-12  Siqi Li  Support more global HYCOM datasets
# 2023-01-13  Siqi Li  Merge the data in type 2
#                      Add the time interval setting
# 2023-01-19  Siqi Li  Set limits based on longitudes and latidues, rather
#                      than index
#==========================================================================


#-----Settings-------------------------------------------------------------
# Starting time (HH:MM yyyy-mm-dd)
time1='00:00 2017-01-01'
# Ending time (HH:MM yyyy-mm-dd)
time2='09:00 2017-01-01'
# Ouput directory
outdir='../hycom'
expt='v57.2'
# Time interval
dt='3 hour'
# Longitude and latitude limits in id (ncks required)
# Use lonlims=() to download data of all the longitude
# Use latlims=() to download data of all the latitude
# Example: lonlims=(-130.88 -127.76); latlims=(47.12 48.84)
lonlims=(-77.9642 -56.8507)
latlims=(31.8400 46.1460)
#--------------------------------------------------------------------------

# Starting time
time=${time1}
t1=`date -u -d "${time}" '+%Y%m%d%H%M'`
t=$t1

# Ending time
t2=`date -u -d "${time2}" '+%Y%m%d%H%M'`

# URL
glb=${expt:0:1}
expt=${expt:1:4}
url0='https://ncss.hycom.org/thredds/ncss/GLB'${glb}'0.08/expt_'${expt}
if [ $expt == '53.X' ]; then
  url0+='/data/YYYY'
fi
url0+='?'
url0+='var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v&'
url0+='disableProjSubset=on&horizStride=1&vertCoord=&addLatLon=true&accept=netcdf4&'

# Deal with longitudes and latitudes
type1="93.0 92.9 92.8"
type2="57.7 57.2 56.3 53.X"
if [[ "$type1" =~ "$expt" ]]; then
  # Set the default lonlims and latlims
  if [ -z "$lonlims" ]; then
    lonlims=(0.0 359.92)
  fi
  if [ -z "$latlims" ]; then
    latlims=(-80.0 90.0)
  fi
  # Convert longitude from [-180 180] to [0 360]
  if (( $(echo "${lonlims[0]} < 0" |bc -l) )); then
    lonlims[0]=`echo "${lonlims[0]} + 360.0" | bc`
    lonlims[1]=`echo "${lonlims[1]} + 360.0" | bc`
  fi
elif [[ "$type2" =~ "$expt" ]]; then
  # Set the default lonlims and latlims
  if [ -z "$lonlims" ]; then
    lonlims=(-180.0 179.92)
  fi
  if [ -z "$latlims" ]; then
    latlims=(-80.0 90.0)
  fi
  # Convert longitude from [0 360] to [-180 180]
  if (( $(echo "${lonlims[1]} > 180.0" |bc -l) )); then
    lonlims[0]=`echo "${lonlims[0]} - 360.0" | bc`
    lonlims[1]=`echo "${lonlims[1]} - 360.0" | bc`
  fi
else
  echo '============================='
  echo '| Unknown dataset: ' $expt '|'
  echo '============================='
  exit
fi

# Loop
while [ $t -le $t2 ]; do

  # Date 
  yyyy=`date -u -d "${time}" '+%Y'`
  mm=`date -u -d "${time}" '+%m'`
  dd=`date -u -d "${time}" '+%d'`
  HH=`date -u -d "${time}" '+%H'`

  echo $yyyy $mm $dd $HH 

  # URL
  url=${url0}
  url+='south='${latlims[0]}'&'
  url+='north='${latlims[1]}'&'
  url+='west='${lonlims[0]}'&'
  url+='east='${lonlims[1]}'&'
  url+='time='${yyyy}'-'${mm}'-'${dd}'T'${HH}'%3A00%3A00Z'
  url=`echo "${url}" | sed -r 's/yyyy/${yyyy}/g'`

  # File name
  fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00.nc

  wget -q ${url} -O ${fout} -o - 
#  echo ${url}
  echo
  echo

  # Advance time
  time=`date -u -d "${time} +${dt}" '+%H:%M %Y-%m-%d'`
  t=`date -u -d "${time}" '+%Y%m%d%H%M '`

done
