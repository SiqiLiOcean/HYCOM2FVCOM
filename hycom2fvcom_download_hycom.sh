#!/bin/bash

#==========================================================================
# Download HYCOM dataset 
#
# https://www.hycom.org/dataserver
#
# Supported datasets
# 1 --- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12° Reanalysis 
#       1994-01-01 to 2015-12-31
# 2 --- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12° Analysis 
#       2014-July to Present
#
# resolution: 0.08 x 0.08 
# time interval: 3 hours
# format: NetCDF
#
# Siqi Li, SMAST
# 2022-11-30
#
# Updates
# 2022-12-16  Siqi Li  Add the settings of longitude/latitude limits
#==========================================================================


#-----Settings-------------------------------------------------------------
# Starting time (HH:MM yyyy-mm-dd)
time1='00:00 2020-01-01'
# Ending time (HH:MM yyyy-mm-dd)
time2='00:00 2020-01-01'
# Ouput directory
outdir='./test'
# Dataset
# 1 --- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12° Reanalysis 
#       1994-01-01 to 2015-12-31
# 2 --- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12° Analysis 
#       2014-July to Present
dataset=2
# Longitude and latitude limits (ncks required)
# Use xlims=() to download data of all the longitude
# Use ylims=() to download data of all the latitude
# Example: xlims=(1470,1865); ylims=(2478,3301)
xlims=(1470,1865)
ylims=(2478,3301)
#--------------------------------------------------------------------------

# Starting time
time=${time1}
t=`date -u -d "${time}" '+%Y%m%d%H%M'`

# Ending time
t2=`date -u -d "${time2}" '+%Y%m%d%H%M'`

# Longitude and Latitude limits
if [ -z "$xlims" ]; then
  xlims_str=''
else
  xlims_str='-d lon,'${xlims}
fi
if [ -z "$ylims" ]; then
  ylims_str=''
else
  ylims_str='-d lat,'${ylims}
fi

# Loop
while [ $t -le $t2 ]; do

  # Date 
  yyyy=`date -u -d "${time}" '+%Y'`
  mm=`date -u -d "${time}" '+%m'`
  dd=`date -u -d "${time}" '+%d'`
  HH=`date -u -d "${time}" '+%H'`
  fyyyy=`date -u -d "${time} -12 hour" '+%Y'`
  fmm=`date -u -d "${time} -12 hour" '+%m'`
  fdd=`date -u -d "${time} -12 hour" '+%d'`
  fHH=`date -u -d "${time} -12 hour" '+%H'`

echo $yyyy $mm $dd $HH --- $fyyyy $fmm $fdd $fHH

  if [ $dataset == '1' ]; then
    # Download the data
    fname=hycom_GLBv0.08_532_${fyyyy}${fmm}${fdd}12_t0${fHH}.nc
    url0=http://data.hycom.org/datasets/GLBv0.08/expt_53.X/data/${fyyyy}/
    url=${url0}/${fname}
    # File names
    fin=${outdir}/${fname}
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00.nc
    if [ -z "$xlims" ] && [ -z "$ylims" ]; then
      wget ${url} -P ${outdir}
      mv ${fin} ${fout}
    else
      ncks -F ${xlims_str} ${ylims_str} ${url} -o ${fout}
      rm -rf $datasets
    fi

  elif [ $dataset == '2' ]; then
    for var in ssh ts3z uv3z; do
      # Download the data
      fname=hycom_glby_930_${fyyyy}${fmm}${fdd}12_t0${fHH}_${var}.nc
      url0=http://data.hycom.org/datasets/GLBy0.08/expt_93.0/data/hindcasts/${fyyyy}/
      url=${url0}/${fname}
      # File names
      fin=${outdir}/${fname}
      fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_${var}.nc
      if [ -z "$xlims" ] && [ -z "$ylims" ]; then
        wget ${url} -P ${outdir} 
        mv ${fin} ${fout}
      else
        ncks -F ${xlims_str} ${ylims_str} ${url} -o ${fout}
	rm -rf $datasets
      fi
    done

  else
    echo 'Wrong dataset type'
    exit
  fi

  # Advance time
  time=`date -u -d "${time} +3 hour" '+%H:%M %Y-%m-%d'`
  t=`date -u -d "${time}" '+%Y%m%d%H%M '`

done
