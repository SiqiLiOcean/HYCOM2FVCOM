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
# Siqi Li, SMAST
# 2022-11-30
#
# Updates
# 2022-12-16  Siqi Li  Add the settings of longitude/latitude limits
# 2023-01-12  Siqi Li  Support more global HYCOM datasets
# 2023-01-13  Siqi Li  Merge the data in type 2
#                      Add the time interval setting
#==========================================================================


#-----Settings-------------------------------------------------------------
# Starting time (HH:MM yyyy-mm-dd)
time1='21:00 2016-12-23'
# Ending time (HH:MM yyyy-mm-dd)
time2='15:00 2016-12-29'
# Ouput directory
outdir='../hycom'
expt='57.2'
#expt='92.8'
# Time interval
dt='3 hour'
# Longitude and latitude limits in id (ncks required)
# Use xlims=() to download data of all the longitude
# Use ylims=() to download data of all the latitude
# Example: xlims=(1470,1865); ylims=(2478,3301)
xlims=(1276,1541)
ylims=(1898,2155)
#xlims=(3526,3791)
#ylims=(1898,2155)
#--------------------------------------------------------------------------

# Starting time
time=${time1}
t1=`date -u -d "${time}" '+%Y%m%d%H%M'`
t=$t1

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

# URL
url0='http://data.hycom.org/datasets/GLBv0.08/expt_'${expt}'/data'

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

  # Check the data type
  if [ $t -eq $t1 ]; then
    type=1
    list=`curl -l -s $url0/$yyyy/`
    for name in $list; do
      if [[ "$name" =~ "_ssh" ]]; then
        type=2
        break
      fi
    done
  fi

echo $yyyy $mm $dd $HH --- $fyyyy $fmm $fdd $fHH

  if [ $type == '1' ]; then
    # Download the zeta, temperature, salinity, u, v
    fname=hycom_GLBv0.08_${expt:0:2}${expt:3:1}_${fyyyy}${fmm}${fdd}12_t0${fHH}.nc
    url=${url0}/${fyyyy}/${fname}
    fin=${outdir}/${fname}
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00.nc
    echo '  --- Download zeta, t, s, u, v'
    wget -q ${url} -P ${outdir} -o - 
    mv ${fin} ${fout}
    echo '  --- File operation'
    if ! [ -z "$xlims" ] || ! [ -z "$ylims" ]; then
      ncks -O -F ${xlims_str} ${ylims_str} ${fout} -o ${fout}
    fi
    ncks -O -x -v tau,water_temp_bottom,salinity_bottom,water_u_bottom,water_v_bottom ${fout} -o ${fout}

  elif [ $type == '2' ]; then

    # Download zeta
    fname=hycom_glbv_${expt:0:2}${expt:3:1}_${fyyyy}${fmm}${fdd}12_t0${fHH}_ssh.nc
    url=${url0}/${fyyyy}/${fname}
    fin=${outdir}/${fname}
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_ssh.nc
    echo '  --- Download zeta'
    wget -q ${url} -P ${outdir} -o - 
    mv ${fin} ${fout}
    echo '  --- File operation'
    if ! [ -z "$xlims" ] || ! [ -z "$ylims" ]; then
      ncks -O -F ${xlims_str} ${ylims_str} ${fout} -o ${fout}
    fi
    ncks -O -x -v tau ${fout} -o ${fout}

    # Download temperature, salinity
    fname=hycom_glbv_${expt:0:2}${expt:3:1}_${fyyyy}${fmm}${fdd}12_t0${fHH}_ts3z.nc
    url=${url0}/${fyyyy}/${fname}
    fin=${outdir}/${fname}
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_ts3z.nc
    echo '  --- Download t, s'
    wget -q ${url} -P ${outdir} -o -
    mv ${fin} ${fout}
    echo '  --- File operation'
    if ! [ -z "$xlims" ] || ! [ -z "$ylims" ]; then
      ncks -O -F ${xlims_str} ${ylims_str} ${fout} -o ${fout}
    fi
    ncks -O -x -v tau,water_temp_bottom,salinity_bottom ${fout} -o ${fout}

    # Download u, v
    fname=hycom_glbv_${expt:0:2}${expt:3:1}_${fyyyy}${fmm}${fdd}12_t0${fHH}_uv3z.nc
    url=${url0}/${fyyyy}/${fname}
    fin=${outdir}/${fname}
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_uv3z.nc
    echo '  --- Download u, v'
    wget -q ${url} -P ${outdir} -o -
    mv ${fin} ${fout}
    if ! [ -z "$xlims" ] || ![ -z "$ylims" ]; then
      ncks -O -F ${xlims_str} ${ylims_str} ${fout} -o ${fout}
    fi  
    echo '  --- File operation'
    ncks -O -x -v tau,water_u_bottom,water_v_bottom ${fout} -o ${fout}

    # Merge the three files into one
    echo '  --- Combine files'
    fssh=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_ssh.nc
    fts=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_ts3z.nc
    fuv=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00_uv3z.nc
    fout=${outdir}/hycom_${yyyy}${mm}${dd}_${HH}00.nc
    ncks -A $fssh $fts
    ncks -A $fts $fuv
    mv $fuv $fout
    rm -rf $fssh
    rm -rf $fts

  fi

  # Advance time
  time=`date -u -d "${time} +${dt}" '+%H:%M %Y-%m-%d'`
  t=`date -u -d "${time}" '+%Y%m%d%H%M '`

done
