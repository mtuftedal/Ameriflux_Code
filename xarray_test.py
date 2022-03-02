import xarray as xr
import os
import numpy as np
import pandas as pd


FLUX_path=(r'Z:/Matt/virtual_machine_shared/Ameriflux/sgpfluxE39.b1/')
FLUX_files=os.listdir(FLUX_path)
for files in sorted(FLUX_files):
    fpath_FLUX=os.path.join(FLUX_path,files)
    ds = xr.open_dataset(fpath_FLUX,engine='netcdf4')
    time=ds['time'].data
    time_bounds=ds['time_bounds_ecor'].data

#-----ECORSF Data. Ameriflux data variable assigned for consistency.------
    
    # Carbon Dioxide (CO2) turbulent flux (no storage correction)
    try:    
        FC=ds['co2_flux'].rename("FC")
    except:
        #FC=np.ones((time.shape))*np.nan
        FC=np.ones((time.shape))*-9999
        FC=xr.DataArray(FC,dims='time')
        FC=FC.rename("FC")
  
    # Methane (CH4) turbulent flux (no storage correction)
    # FCH4 Not in our data
    
    # Carbon Dioxide (CO2) mole fraction in wet air
    try:    
        CO2=ds['co2_molar_fraction'].rename("CO2")
    except:
        #CO2=np.ones((time.shape))*np.nan
        CO2=np.ones((time.shape))*-9999
        CO2=xr.DataArray(CO2,dims='time')
        CO2=CO2.rename("CO2")

    # Carbon Dioxide (CO2) in mole fraction of dry air
    try:    
        CO2_MIXING_RATIO=ds['co2_mixing_ratio'].rename("CO2_MIXING_RATIO")

    except:
        #CO2_MIXING_RATIO=np.ones((time.shape))*np.nan
        CO2_MIXING_RATIO=np.ones((time.shape))*-9999
        CO2_MIXING_RATIO=xr.DataArray(CO2_MIXING_RATIO,dims='time')
        CO2_MIXING_RATIO=CO2_MIXING_RATIO.rename("CO2_MIXING_RATIO")
        
    # Water (H2O) vapor in mole fraction of wet air
    try:    
        H2O=ds['h2o_molar_fraction'].rename("H20")
    except:
        #H2O=np.ones((time.shape))*np.nan
        H2O=np.ones((time.shape))*-9999
        H2O=xr.DataArray(H2O,dims='time')
        H2O=H2O.rename("H2O")    
        
    # Water (H2O) vapor in mole fraction of dry air
    try:    
        H2O_MIXING_RATIO=ds['h2o_mixing_ratio'].rename("H2O_MIXING_RATIO")
    except:
        #H2O_MIXING_RATIO=np.ones((time.shape))*np.nan
        H2O_MIXING_RATIO=np.ones((time.shape))*-9999
        H2O_MIXING_RATIO=xr.DataArray(H2O_MIXING_RATIO,dims='time')
        H2O_MIXING_RATIO=H2O_MIXING_RATIO.rename("H2O_MIXING_RATIO") 
    
    # Methane (CH4) mole fraction in wet air
    try:    
        CH4=ds['ch4_molar_fraction'].rename("CH4")
    except:
        #CH4=np.ones((time.shape))*np.nan
        CH4=np.ones((time.shape))*-9999
        CH4=xr.DataArray(CH4,dims='time')
        CH4=CH4.rename("CH4") 
          
    # Methane (CH4) in mole fraction of dry air    
    try:    
        CH4_MIXING_RATIO=ds['ch4_mixing_ratio'].rename("CH4_MIXING_RATIO")
    except:
        #CH4_MIXING_RATIO=np.ones((time.shape))*np.nan
        CH4_MIXING_RATIO=np.ones((time.shape))*-9999
        CH4_MIXING_RATIO=xr.DataArray(CH4_MIXING_RATIO,dims='time')
        CH4_MIXING_RATIO=CH4_MIXING_RATIO.rename("CH4_MIXING_RATIO")  
        
    # Momentum flux
    try:    
        TAU=ds['momentum_flux'].rename("TAU")
    except:
        #TAU=np.ones((time.shape))*np.nan
        TAU=np.ones((time.shape))*-9999
        TAU=xr.DataArray(TAU,dims='time')
        TAU=TAU.rename("TAU") 
 
    # Sensible heat turbulent flux (no storage correction)
    try:
        H=ds['sensible_heat_flux'].rename("H")
    except:
        #H=np.ones((time.shape))*np.nan
        H=np.ones((time.shape))*-9999
        H=xr.DataArray(H,dims='time')
        H=H.rename("H") 

    # Latent heat turbulent flux (no storage correction)
    try:
        LE=ds['latent_flux'].rename("LE")
    except:
        #LE=np.ones((time.shape))*np.nan
        LE=np.ones((time.shape))*-9999
        LE=xr.DataArray(LE,dims='time')
        LE=LE.rename("LE")         

    # Air temperature given in K then converted to degC
    try:
        TA=ds['air_temperature'].rename("TA")
        TA=TA-273.15
    except:
        #TA=np.ones((time.shape))*np.nan
        TA=np.ones((time.shape))*-9999
        TA=xr.DataArray(TA,dims='time')
        TA=TA.rename("TA") 
        
    # Atmospheric pressure
    try:
        PA=ds['air_pressure'].rename("PA")
    except:
        #PA=np.ones((time.shape))*np.nan
        PA=np.ones((time.shape))*-9999
        PA=xr.DataArray(PA,dims='time')
        PA=PA.rename("PA") 
       
    # Relative humidity, range 0-100
    try:
        RH=ds['relative_humidity'].rename("RH")
    except:
        #RH=np.ones((time.shape))*np.nan
        RH=np.ones((time.shape))*-9999
        RH=xr.DataArray(RH,dims='time')
        RH=RH.rename("RH") 
         
    # Sonic temperature
    try:
        T_SONIC=ds['sonic_temperature'].rename("T_SONIC")
        T_SONIC=T_SONIC-273.15
    except:
        #T_SONIC=np.ones((time.shape))*np.nan
        T_SONIC=np.ones((time.shape))*-9999
        T_SONIC=xr.DataArray(T_SONIC,dims='time')
        T_SONIC=T_SONIC.rename("T_SONIC") 
      
    # Vapor Pressure Deficit and is converted to hPa from kPa
    try:
        VPD=ds['water_vapor_pressure_deficit'].rename("VPD")
        VPD=VPD*10
    except:
        #VPD=np.ones((time.shape))*np.nan
        VPD=np.ones((time.shape))*-9999
        VPD=xr.DataArray(VPD,dims='time')
        VPD=VPD.rename("VPD") 

    # Monin-Obukhov length
    try:
        MO_LENGTH=ds['Monin_Obukhov_length'].rename("MO_LENGTH")
    except:
        #MO_LENGTH=np.ones((time.shape))*np.nan
        MO_LENGTH=np.ones((time.shape))*-9999 
        MO_LENGTH=xr.DataArray(MO_LENGTH,dims='time')
        MO_LENGTH=MO_LENGTH.rename("MO_LEGNTH") 

    # Monin-Obukhov Stability parameter
    try:
        ZL=ds['Monin_Obukhov_stability_parameter'].rename("ZL")
    except:
        #ZL=np.ones((time.shape))*np.nan
        ZL=np.ones((time.shape))*-9999         
        ZL=xr.DataArray(ZL,dims='time')
        ZL=ZL.rename("ZL") 
        
    # Wind speed
    try:
        WS=ds['mean_wind'].rename("WS")
    except:
        #WS=np.ones((time.shape))*np.nan
        WS=np.ones((time.shape))*-9999
        WS=xr.DataArray(WS,dims='time')
        WS=WS.rename("WS") 
    
    # Wind direction
    try:
        WD=ds['wind_direction_from_north'].rename("WD")
    except:
        #WD=np.ones((time.shape))*np.nan
        WD=np.ones((time.shape))*-9999
        WD=xr.DataArray(WD,dims='time')
        WD=WD.rename("WD") 
    
    # Friction velocity
    try:
        USTAR=ds['friction_velocity'].rename("USTAR")
    except:
        #USTAR=np.ones((time.shape))*np.nan
        USTAR=np.ones((time.shape))*-9999
        USTAR=xr.DataArray(USTAR,dims='time')
        USTAR=USTAR.rename("USTAR") 
        
    # Maximum WS in the averaging period
    try:
        WS_MAX=ds['maximum_instantaneous_wind_speed'].rename("WS_MAX")
    except:
        #WS_MAX=np.ones((time.shape))*np.nan
        WS_MAX=np.ones((time.shape))*-9999
        WS_MAX=xr.DataArray(WS_MAX,dims='time')
        WS_MAX=WS_MAX.rename("WS_MAX") 


#-----SEBS Data. Ameriflux data variable assigned for consistency.------
    # Shortwave radiation, incoming
    try:
        SW_IN=ds['down_short_hemisp'].rename('SW_IN')
    except:
        #SW_IN=np.ones((time.shape))*np.nan
        SW_IN=np.ones((time.shape))*-9999
        SW_IN=xr.DataArray(SW_IN,dims='time')
        SW_IN=SW_IN.rename("SW_IN") 

    # Shortwave radiation, outgoing
    try:
        SW_OUT=ds['up_short_hemisp'].rename("SW_OUT")
    except:
        #SW_OUT=np.ones((time.shape))*np.nan
        SW_OUT=np.ones((time.shape))*-9999
        SW_OUT=xr.DataArray(SW_OUT,dims='time')
        SW_OUT=SW_OUT.rename("SW_OUT") 
        
    # Longwave radiation, incoming
    try:
        LW_IN=ds['down_long'].rename("LW_IN")
    except:
        #LW_IN=np.ones((time.shape))*np.nan
        LW_IN=np.ones((time.shape))*-9999        
        LW_IN=xr.DataArray(LW_IN,dims='time')
        LW_IN=H2O.rename("LW_IN") 
        
    # Longwave radiation, outgoing
    try:
        LW_OUT=ds['up_long'].rename("LW_OUT")
    except:
        #LW_OUT=np.ones((time.shape))*np.nan
        LW_OUT=np.ones((time.shape))*-9999
        LW_OUT=xr.DataArray(LW_OUT,dims='time')
        LW_OUT=LW_OUT.rename("LW_OUT") 
        
    # Albedo, range 0-100
    try:
        ALB=ds['albedo'].rename("ALB")
    except:
        #ALB=np.ones((time.shape))*np.nan
        ALB=np.ones((time.shape))*-9999
        ALB=xr.DataArray(ALB,dims='time')
        ALB=ALB.rename("ALB") 

    # Net radiation
    try:
        NETRAD=ds['net_radiation'].rename("NETRAD")
    except:
        #NETRAD=np.ones((time.shape))*np.nan
        NETRAD=np.ones((time.shape))*-9999
        NETRAD=xr.DataArray(NETRAD,dims='time')
        NETRAD=NETRAD.rename("NETRAD") 
        
    # Soil heat flux
    try:
        G_1_1_1=ds['surface_soil_heat_flux_1'].rename("G_1_1_1")
    except:
        #G_1_1_1=np.ones((time.shape))*np.nan
        G_1_1_1=np.ones((time.shape))*-9999
        G_1_1_=xr.DataArray(G_1_1_1,dims='time')
        G_1_1_1=G_1_1_1.rename("G_1_1_1") 
        
    # Soil heat flux
    try:
        G_1_1_2=ds['surface_soil_heat_flux_2'].rename("G_1_1_2")
    except:
        #G_1_1_2=np.ones((time.shape))*np.nan
        G_1_1_2=np.ones((time.shape))*-9999
        G_1_1_2=xr.DataArray(G_1_1_2,dims='time')
        G_1_1_2=G_1_1_2.rename("G_1_1_2") 
        
     # Soil heat flux
    try:
        G_1_1_3=ds['surface_soil_heat_flux_3'].rename("G_1_1_3")
    except:
        #G_1_1_3=np.ones((time.shape))*np.nan
        G_1_1_3=np.ones((time.shape))*-9999
        G_1_1_3=xr.DataArray(G_1_1_3,dims='time')
        G_1_1_3=G_1_1_3.rename("G_1_1_3") 

    # Soil heat flux average. Variable is tested to ensure no missing value
    # codes skew the averaging.Two data points are pulled to make sure data
    # exists in the file. 
    count=0
    if G_1_1_1[0] != -9999 and G_1_1_1[30] !=-9999:
        G1=G_1_1_1
        count=count+1
    else:
        G1=G_1_1_1*0
        
    if G_1_1_2[0] != -9999 and G_1_1_2[30] !=-9999:
        G2=G_1_1_2
        count=count+1
    else:
        G2=G_1_1_2*0
        
    if G_1_1_3[0] != -9999 and G_1_1_3[30] !=-9999:
        G3=G_1_1_3
        count=count+1
    else:
        G3=G_1_1_3*0
    
    G_1_1_A=(G1 + G2 + G2)/count
    G_1_1_A=G_1_1_A.rename("G_1_1_A") 
    
    # Soil temperature
    try:
        TS_1_1_1=ds['soil_temp_1'].rename("TS_1_1_1")
    except:
        #TS_1_1_1=np.ones((time.shape))*np.nan
        TS_1_1_1=np.ones((time.shape))*-9999
        TS_1_1_1=xr.DataArray(TS_1_1_1,dims='time')
        TS_1_1_3=TS_1_1_1.rename("TS_1_1_1") 

    # Soil temperature
    try:
        TS_1_1_2=ds['soil_temp_2'].rename("TS_1_1_2")
    except:
        #TS_1_1_2=np.ones((time.shape))*np.nan
        TS_1_1_2=np.ones((time.shape))*-9999
        TS_1_1_2=xr.DataArray(TS_1_1_2,dims='time')
        TS_1_1_2=TS_1_1_2.rename("TS_1_1_2") 
    
    # Soil temperature
    try:
        TS_1_1_3=ds['soil_temp_3'].rename("TS_1_1_3")
    except:
        #TS_1_1_3=np.ones((time.shape))*np.nan
        TS_1_1_3=np.ones((time.shape))*-9999
        TS_1_1_3=xr.DataArray(TS_1_1_3,dims='time')
        TS_1_1_3=TS_1_1_3.rename("TS_1_1_3") 

    # Soil temperature average. Variables are tested to verify they exist
    # and ensure no missing value codes skew the data. 
    count=0
    if TS_1_1_1[0] != -9999 and TS_1_1_1[30] !=-9999:
        TS1_1=TS_1_1_1
        count=count+1
    else:
        TS1=TS_1_1_1*0
    if TS_1_1_2[0] != -9999 and TS_1_1_2[30] !=-9999:
        TS1_2=TS_1_1_2
        count=count+1
    else:
        TS1_2=TS_1_1_2*0
    if TS_1_1_3[0] != -9999 and TS_1_1_3[30] !=-9999:
        TS1_3=TS_1_1_3
        count=count+1
    else:
        TS1_3=TS_1_1_3*0
    
    TS_1_1_A=(TS1_1 + TS1_2 + TS1_3)/count
    TS_1_1_A=TS_1_1_A.rename("TS_1_1_A") 


#-----AMC Data. Ameriflux data variable assigned for consistency.------    
    # Soil Temperature from AMC probe depth -36.8 cm
    try:
        TS_2_2_1=ds['temp_1'].rename("TS_2_2_1")
    except:
        #TS_2_2_1=np.ones((time.shape))*np.nan
        TS_2_2_1=np.ones((time.shape))*-9999 
        TS_2_2_1=xr.DataArray(TS_2_2_1,dims='time')
        TS_2_2_1=TS_2_2_1.rename("TS_2_2_1") 

    # Soil temperature from AMC probe depth -14 cm
    try:
        TS_2_1_1=ds['temp_2'].rename("TS_2_1_1")
    except:
        #TS_2_1_1=np.ones((time.shape))*np.nan
        TS_2_1_1=np.ones((time.shape))*-9999
        TS_2_1_1=xr.DataArray(TS_2_1_1,dims='time')
        TS_2_1_1=TS_2_1_1.rename("TS_2_1_1")          

    # Soil temperature from AMC probe depth -35.6 cm
    try:
        TS_2_2_2=ds['temp_3'].rename("TS_2_2_2")
    except:
        #TS_2_2_2=np.ones((time.shape))*np.nan
        TS_2_2_2=np.ones((time.shape))*-9999
        TS_2_2_2=xr.DataArray(TS_2_2_2,dims='time')
        TS_2_2_2=TS_2_2_2.rename("TS_2_2_2") 
        
    # Soil temperature from AMC probe depth -14 cm
    try:
        TS_2_1_2=ds['temp_4'].rename("TS_2_1_2")
    except:
        #TS_2_1_2=np.ones((time.shape))*np.nan
        TS_2_1_2=np.ones((time.shape))*-9999
        TS_2_1_2=xr.DataArray(TS_2_1_2,dims='time')
        TS_2_1_2=TS_2_1_2.rename("TS_2_1_2") 

    # Soil temperature from AMC probe depth -35.6 cm
    try:
        TS_2_2_3=ds['temp_5'].rename("TS_2_2_3")
    except:
        #TS_2_2_3=np.ones((time.shape))*np.nan
        TS_2_2_3=np.ones((time.shape))*-9999
        TS_2_2_3=xr.DataArray(TS_2_2_3,dims='time')
        TS_2_2_3=TS_2_2_3.rename("TS_2_2_3") 

    # Soil temperature from AMC probe depth -14 cm
    try:
        TS_2_1_3=ds['temp_6'].rename("TS_2_1_3")
    except:
        #TS_2_1_3=np.ones((time.shape))*np.nan
        TS_2_1_3=np.ones((time.shape))*-9999
        TS_2_1_3=xr.DataArray(TS_2_1_3,dims='time')
        TS_2_1_3=TS_2_1_3.rename("TS_2_1_3") 

    # Soil temperature from AMC probe depth -34.3 cm
    try:
        TS_2_2_4=ds['temp_7'].rename("TS_2_2_4")
    except:
        #TS_2_2_4=np.ones((time.shape))*np.nan
        TS_2_2_4=np.ones((time.shape))*-9999 
        TS_2_2_4=xr.DataArray(TS_2_2_4,dims='time')
        TS_2_2_4=TS_2_2_4.rename("TS_2_2_4") 
        
    # Soil temperature from AMC probe depth -15.2 cm
    try:
        TS_2_1_4=ds['temp_8'].rename("TS_2_1_4")
    except:
        #TS_2_1_4=np.ones((time.shape))*np.nan
        TS_2_1_4=np.ones((time.shape))*-9999
        TS_2_1_4=xr.DataArray(TS_2_1_4,dims='time')
        TS_2_1_4=TS_2_1_4.rename("TS_2_1_4") 
        
     # Soil temperature from AMC probe depth -34.9
    try:
        TS_2_2_5=ds['temp_9'].rename("TS_2_2_5")
    except:
        #TS_2_2_5=np.ones((time.shape))*np.nan
        TS_2_2_5=np.ones((time.shape))*-9999
        TS_2_2_5=xr.DataArray(TS_2_2_5,dims='time')
        TS_2_2_5=TS_2_2_5.rename("TS_2_2_5")
     
    # Soil temperature from AMC probe depth -14 cm
    try:
        TS_2_1_5=ds['temp_10'].rename("TS_2_1_5")
    except:
        #TS_2_1_5=np.ones((time.shape))*np.nan
        TS_2_1_5=np.ones((time.shape))*-9999
        TS_2_1_5=xr.DataArray(TS_2_1_5,dims='time')
        TS_2_1_5=TS_2_1_5.rename("TS_2_1_5")
        
    # Soil temperature from AMC probe depth -34.3 cm
    try:
        TS_2_2_6=ds['temp_11'].rename("TS_2_2_6")
    except:
        #TS_2_2_6=np.ones((time.shape))*np.nan
        TS_2_2_6=np.ones((time.shape))*-9999
        TS_2_2_6=xr.DataArray(TS_2_2_6,dims='time')
        TS_2_2_6=TS_2_2_6.rename("TS_2_2_6")
   
    # Soil temperature from AMC probe depth -16.5 cm
    try:
        TS_2_1_6=ds['temp_12'].rename("TS_2_1_6")
    except:
        #TS_2_1_6=np.ones((time.shape))*np.nan
        TS_2_1_6=np.ones((time.shape))*-9999
        TS_2_1_6=xr.DataArray(TS_2_1_6,dims='time')
        TS_2_1_6=TS_2_1_6.rename("TS_2_1_4")
        
    # Temperature averaging for probes closer (-14 to -16 cm) to ground surface.
    # The code checks to verify the data exists and for the averaging. If the 
    # data is -9999, the value is multipled by zero. 
    count=0
    if TS_2_1_1[0] != -9999 and TS_2_1_1[30] !=-9999:
        TS_2_1=TS_2_1_1
        count=count+1
    else:
        TS_2_1=TS_2_1_1*0
        
    if TS_2_1_2[0] != -9999 and TS_2_1_2[30] !=-9999:
        TS_2_2=TS_2_1_2
        count=count+1
    else:
        TS_2_2=TS_2_1_2*0
        
    if TS_2_1_3[0] != -9999 and TS_2_1_3[30] !=-9999:
        TS_2_3=TS_2_1_3
        count=count+1
    else:
        TS_2_3=TS_2_1_3*0
        
    if TS_2_1_4[0] != -9999 and TS_2_1_4[30] !=-9999:
        TS_2_4=TS_2_1_4
        count=count+1
    else:
        TS_2_4=TS_2_1_4*0
        
    if TS_2_1_5[0] != -9999 and TS_2_1_5[30] !=-9999:
        TS_2_5=TS_2_1_5
        count=count+1
    else:
        TS_2_5=TS_2_1_5*0
        
    if TS_2_1_6[0] != -9999 and TS_2_1_6[30] !=-9999:
        TS_2_6=TS_2_1_6
        count=count+1
    else:
        TS_2_6=TS_2_1_6*0
    
    # The averaging is done here and then renames the variable in xarray. 
    TS_2_1_A=(TS_2_1 + TS_2_2 + TS_2_3 + TS_2_4 + TS_2_5 + TS_2_6)/count
    TS_2_1_A=TS_2_1_A.rename("TS_2_1_A")
    
    # Temperature averaging for probes below -30 cm. Same averaging method
    # is applied here as previous steps. 
    count=0
    if TS_2_2_1[0] != -9999 and TS_2_2_1[30] !=-9999:
        TS_3_1=TS_2_2_1
        count=count+1
    else:
        TS_3_1=TS_2_2_1*0
        
    if TS_2_2_2[0] != -9999 and TS_2_2_2[30] !=-9999:
        TS_3_2=TS_2_2_2
        count=count+1
    else:
        TS_3_2=TS_2_2_2*0
        
    if TS_2_2_3[0] != -9999 and TS_2_2_3[30] !=-9999:
        TS_3_3=TS_2_2_3
        count=count+1
    else:
        TS_3_3=TS_2_2_3*0
        
    if TS_2_2_4[0] != -9999 and TS_2_2_4[30] !=-9999:
        TS_3_4=TS_2_2_4
        count=count+1
    else:
        TS_3_4=TS_2_2_4*0
        
    if TS_2_2_5[0] != -9999 and TS_2_2_5[30] !=-9999:
        TS_3_5=TS_2_2_5
        count=count+1
    else:
        TS_3_5=TS_2_2_5*0
        
    if TS_2_2_6[0] != -9999 and TS_2_2_6[30] !=-9999:
        TS_3_6=TS_2_2_6
        count=count+1
    else:
        TS_3_6=TS_2_2_6*0
    
    # Averaging and renaming the variable. 
    TS_2_2_A=(TS_3_1 + TS_3_2 + TS_3_3 + TS_3_4 + TS_3_5 + TS_3_6)/count   
    TS_2_2_A=TS_2_2_A.rename("TS_2_2_A")
    
    # Photosynthetic photon flux density, incoming
    try:
        PPFD_IN=ds['par_inc'].rename("PPFD_IN")
    except:
        #PPFD_IN=np.ones((time.shape))*np.nan
        PPFD_IN=np.ones((time.shape))*-9999
        PPFD_IN=xr.DataArray(PPFD_IN,dims='time')
        PPFD_IN=PPFD_IN.rename("PPFD_IN")
        
    # Photosynthetic photon flux density, outgoing
    try:
        PPFD_OUT=ds['par_ref'].rename("PPFD_OUT")
    except:
        #PPFD_OUT=np.ones((time.shape))*np.nan
        PPFD_OUT=np.ones((time.shape))*-9999
        PPFD_OUT=xr.DataArray(PPFD_OUT,dims='time')
        PPFD_OUT=PPFD_OUT.rename("PPFD_OUT")
    
    # Soil water content (volumetric), range 0-100
    # Soil water content at -36.8 cm
    try:
        SWC_1_2_1=ds['vwc_1'].rename("SWC_1_2_1")
    except:
        #SWC_1_2_1=np.ones((time.shape))*np.nan
        SWC_1_2_1=np.ones((time.shape))*-9999 
        SWC_1_2_1=xr.DataArray(SWC_1_2_1,dims='time')
        SWC_1_2_1=SWC_1_2_1.rename("SWC_1_2_1")
        
    # Soil water content at -14 cm
    try:
        SWC_1_1_1=ds['vwc_2'].rename("SWC_1_1_1")
    except:
        #SWC_1_1_1=np.ones((time.shape))*np.nan
        SWC_1_1_1=np.ones((time.shape))*-9999
        SWC_1_1_1=xr.DataArray(SWC_1_1_1,dims='time')
        SWC_1_1_1=SWC_1_1_1.rename("SWC_1_1_1")
    
    # Soil water content at -35.6 cm
    try:
        SWC_1_2_2=ds['vwc_3'].rename("SWC_1_2_2")
    except:
        #SWC_1_2_2=np.ones((time.shape))*np.nan
        SWC_1_2_2=np.ones((time.shape))*-9999
        SWC_1_2_2=xr.DataArray(SWC_1_2_2,dims='time')
        SWC_1_2_2=SWC_1_2_2.rename("SWC_1_2_2")
        
    # Soil water content at -14 cm
    try:
        SWC_1_1_2=ds['vwc_4'].rename("SWC_1_1_2")
    except:
        #SWC_1_1_2=np.ones((time.shape))*np.nan
        SWC_1_1_2=np.ones((time.shape))*-9999
        SWC_1_1_2=xr.DataArray(SWC_1_1_2,dims='time')
        SWC_1_1_2=SWC_1_1_2.rename("SWC_1_1_2")
        
    # Soil water content at -35.6 cm    
    try:
        SWC_1_2_3=ds['vwc_5'].rename("SWC_1_2_3")
    except:
        #SWC_1_2_3=np.ones((time.shape))*np.nan
        SWC_1_2_3=np.ones((time.shape))*-9999
        SWC_1_2_3=xr.DataArray(SWC_1_2_3,dims='time')
        SWC_1_2_3=SWC_1_2_3.rename("SWC_1_2_3")
        
    # Soil water content at -14 cm    
    try:
        SWC_1_1_3=ds['vwc_6'].rename("SWC_1_1_3")
    except:
        #SWC_1_1_3=np.ones((time.shape))*np.nan
        SWC_1_1_3=np.ones((time.shape))*-9999
        SWC_1_1_3=xr.DataArray(SWC_1_1_3,dims='time')
        SWC_1_1_3=SWC_1_1_3.rename("SWC_1_1_3")
            
    # Soil water content at -34.3 cm
    try:
        SWC_1_2_4=ds['vwc_7'].rename("SWC_1_2_4")
    except:
        #SWC_1_2_4=np.ones((time.shape))*np.nan
        SWC_1_2_4=np.ones((time.shape))*-9999
        SWC_1_2_4=xr.DataArray(SWC_1_2_4,dims='time')
        SWC_1_2_4=SWC_1_2_4.rename("SWC_1_2_4")
        
    # Soil water content at -15.2 cm    
    try:
        SWC_1_1_4=ds['vwc_8'].rename("SWC_1_1_4")
    except:
        #SWC_1_1_4=np.ones((time.shape))*np.nan
        SWC_1_1_4=np.ones((time.shape))*-9999
        SWC_1_1_4=xr.DataArray(SWC_1_1_4,dims='time')
        SWC_1_1_4=SWC_1_1_4.rename("SWC_1_1_4")
                    
    # Soil water content at -34.9 cm
    try:
        SWC_1_2_5=ds['vwc_9'].rename("SWC_1_2_5")
    except:
        #SWC_1_2_5=np.ones((time.shape))*np.nan
        SWC_1_2_5=np.ones((time.shape))*-9999
        SWC_1_2_5=xr.DataArray(SWC_1_2_5,dims='time')
        SWC_1_2_5=SWC_1_2_5.rename("SWC_1_2_5")

    # Soil water content at -14 cm    
    try:
        SWC_1_1_5=ds['vwc_10'].rename("SWC_1_1_5")
    except:
        #SWC_1_1_5=np.ones((time.shape))*np.nan
        SWC_1_1_5=np.ones((time.shape))*-9999
        SWC_1_1_5=xr.DataArray(SWC_1_1_5,dims='time')
        SWC_1_1_5=SWC_1_1_5.rename("SWC_1_1_5")

    # Soil water content at -34.3 cm
    try:
        SWC_1_2_6=ds['vwc_11'].rename("SWC_1_2_6")
    except:
        #SWC_1_2_6=np.ones((time.shape))*np.nan
        SWC_1_2_6=np.ones((time.shape))*-9999
        SWC_1_2_6=xr.DataArray(SWC_1_2_6,dims='time')
        SWC_1_2_6=SWC_1_2_6.rename("SWC_1_2_6")

    # Soil water content at -16.5 cm    
    try:
        SWC_1_1_6=ds['vwc_12'].rename("SWC_1_1_6")
    except:
        #SWC_1_1_6=np.ones((time.shape))*np.nan
        SWC_1_1_6=np.ones((time.shape))*-9999
        SWC_1_1_6=xr.DataArray(SWC_1_1_6,dims='time')
        SWC_1_1_6=SWC_1_1_6.rename("SWC_1_1_6")

    # Soil water content averaging for probes closer (-14 to -16 cm) to ground surface. 
    count=0
    if SWC_1_1_1[0] != -9999 and SWC_1_1_1[30] !=-9999:
        SWC_1_1=SWC_1_1_1
        count=count+1
    else:
        SWC_1_1=SWC_1_1_1*0
        
    if SWC_1_1_2[0] != -9999 and SWC_1_1_2[30] !=-9999:
        SWC_1_2=SWC_1_1_2
        count=count+1
    else:
        SWC_1_2=SWC_1_1_2*0
        
    if SWC_1_1_3[0] != -9999 and SWC_1_1_3[30] !=-9999:
        SWC_1_3=SWC_1_1_3
        count=count+1
    else:
        SWC_1_3=SWC_1_1_3*0
        
    if SWC_1_1_4[0] != -9999 and SWC_1_1_4[30] !=-9999:
        SWC_1_4=SWC_1_1_4
        count=count+1
    else:
        SWC_1_4=SWC_1_1_4*0
        
    if SWC_1_1_5[0] != -9999 and SWC_1_1_5[30] !=-9999:
        SWC_1_5=SWC_1_1_5
        count=count+1
    else:
        SWC_1_5=SWC_1_1_5*0
        
    if SWC_1_1_6[0] != -9999 and SWC_1_1_6[30] !=-9999:
        SWC_1_6=SWC_1_1_6
        count=count+1
    else:
        SWC_1_6=SWC_1_1_6*0
        
    SWC_1_1_A=(SWC_1_1 + SWC_1_2 + SWC_1_3 + SWC_1_4 + SWC_1_5 + SWC_1_6)/count 
    SWC_1_1_A=SWC_1_1_A.rename("SWC_1_1_A")
    
    # Soil water content averaging for probes between -34.0 cm and -36.0 cm . 
    count=0
    if SWC_1_2_1[0] != -9999 and SWC_1_2_1[30] !=-9999:
        SWC_2_1=SWC_1_2_1
        count=count+1
    else:
        SWC_2_1=SWC_1_2_1*0
        
    if SWC_1_2_2[0] != -9999 and SWC_1_2_2[30] !=-9999:
        SWC_2_2=SWC_1_2_2
        count=count+1
    else:
        SWC_2_2=SWC_1_2_2*0
        
    if SWC_1_2_3[0] != -9999 and SWC_1_2_3[30] !=-9999:
        SWC_2_3=SWC_1_2_3
        count=count+1
    else:
        SWC_2_3=SWC_1_2_3*0
        
    if SWC_1_2_4[0] != -9999 and SWC_1_2_4[30] !=-9999:
        SWC_2_4=SWC_1_2_4
        count=count+1
    else:
        SWC_2_4=SWC_1_2_4*0
        
    if SWC_1_2_5[0] != -9999 and SWC_1_2_5[30] !=-9999:
        SWC_2_5=SWC_1_2_5
        count=count+1
    else:
        SWC_2_5=SWC_1_2_5*0
        
    if SWC_1_2_6[0] != -9999 and SWC_1_2_6[30] !=-9999:
        SWC_2_6=SWC_1_2_6
        count=count+1
    else:
        SWC_2_6=SWC_1_2_6*0
        
    SWC_1_2_A=(SWC_2_1 + SWC_2_2 + SWC_2_3 + SWC_2_4 + SWC_2_5 + SWC_2_6)/count
    SWC_1_2_A=SWC_1_2_A.rename("SWC_1_2_A")
       



