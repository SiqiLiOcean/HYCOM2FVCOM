# HYCOM2FVCOM

Siqi Li, SMAST
2022-12-06

### Use HYCOM dataset to run FVCOM case
    Create initial TS file
    Create nesting forcing file

### The package including:
    hycom2fvcom_download_hycom.sh
    hycom2fvcom_iniTS_create.m
    hycom2fvcom_nesting_select.m
    hycom2fvcom_nesting_tide_TMD.m
    hycom2fvcom_nesting_tide_user.m
    hycom2fvcom_nesting_create.m


### Required MATLAB toolbox
+ matFVCOM
https://github.com/SiqiLiOcean/matFVCOM
+ matFigure
https://github.com/SiqiLiOcean/matFigure
+ t_tide
https://www.eoas.ubc.ca/~rich/t_tide/t_tide_v1.4beta.zip
+ TMD
https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5
    
### Required data
Ocean tide model data is required to extract the tide constituents. 
+ Visit https://www.tpxo.net/global, and register there following 'HOW TO REGISTER'.
+ Write the 'load' file, containing three lines for the paths and names of the downloaded 'h_', 'u_', and 'grid_'. For example,
            ~/tools/TMD/DATA/h_tpxo9.v5a  
            ~/tools/TMD/DATA/u_tpxo9.v5a  
            ~/tools/TMD/DATA/grid_tpxo9.v5a
+ Set the 'load' file path and name to 'Model' in hycom2fvcom_nesting_extract_tide.m.


### Codes
-------------------------- Initial TS file --------------------------------
1. Create the intial TS file (hycom2fvcom_iniTS_create.m)


-------------------------- Nesting forcing --------------------------------
1. Download HYCOM data (hycom2fvcom_download_hycom.sh)
    bash hycom2fvcom_download_hycom.sh
      OR
   download data manually:
    + Select the proper HYCOM dataset from https://www.hycom.org/dataserver,
      based on the region, resolution, and time period
    + Click the selected dataset, and then click 'ACCESS DATA HERE' on the
      top.
    + Click the right directory. For example:
      + for Japan region of 2002
       --> GLBv0.08/expt_53.X with Sea Ice (Jan-01-1994 to Dec-30-2015)/
           http://data.hycom.org/datasets/GLBv0.08/expt_53.X/data/2002/   
           + hycom_GLBv0.08_532_2002010112_t000.nc
           + hycom_GLBv0.08_532_2002010112_t003.nc
           + ...
      + for Japan region of 2020
       --> GLBy0.08/expt_93.0 (Dec-04-2018 to Present + FORECASTS + ice + sur)/
           http://data.hycom.org/datasets/GLBy0.08/expt_93.0/data/hindcasts/2020/   
           + hycom_glby_930_2020010112_t000_ssh.nc
           + hycom_glby_930_2020010112_t000_ts3z.nc
           + hycom_glby_930_2020010112_t000_uv3z.nc
           + hycom_glby_930_2020010112_t003_ssh.nc
           + hycom_glby_930_2020010112_t003_ts3z.nc
           + hycom_glby_930_2020010112_t003_uv3z.nc
           + ...

2. Set the nesting layers (hycom2fvcom_nesting_select.m)
    Select nesting nodes and cells based on the grid, obc, and the nesting layers

3. Extract the tide harmonic components.
   There are two ways to create the tide components for the FVCOM nesting. One is to use TPXO dataset via TMD toolbox. The other is to use user-defined tide components.
   + a) Use TPXO dataset via TMD toolbox (hycom2fvcom_nesting_tide_TMD.m).
     + download the tide model from
       https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
     + extract the tide harmonic coefficients on the nesting nodes and cells
   + b) User defined (hycom2fvcom_nesting_tide_user.m)
   
4. Write out the nesting forcing file (hycom2fvcom_nesting_create.m)

