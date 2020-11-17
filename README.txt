**** BND_sim ****

Description:
This is a MATLAB code that estimates bubble nucleation during Plinian silicic eruptions. 

The Input variable is a table that contains the following parameter:
1- Water saturation pressure in MPa. The acceptable range is 50-300
2- Temperature in C. The acceptable range is 700-1100.
3- Crystal content for viscosity. The acceptable range is 0-1.
4- Heterogeneous nucleation contact angle in degrees. The acceptable range is 0-180.
5- Mass discharge rate in Kg/s.
6- Conduit radius. 

The output file is a structure file that contains:
1- Time in second
2- Melt pressure in Pa.
3- Water concentration.
4- Gas volume fraction
5- Nucleation rate m^{-3}s^{-1}
5- Bubble number density in m^{-3}
6- Decompression rate in MPa/s


The code uses interpolated variables stored in "magma_properties.mat" to calculate the following properties of magma:
1- Diffusivity as a function of water concentration, pressure, and temperature (Zhang, Y., & Behrens, H. (2000), Chemical Geology).
2- Water fugacity as a function of saturation pressure and temperature (Flowers, G. C. (1979) Contributions to Mineralogy and Petrology).
3- Solubility of water in rhyolite as a function of pressure and temperature (Liu, Y., Zhang, Y., & Behrens, H. (2005), Geothermal Research).
4- Viscosity of rhyolite temperature and water concentration (Hui, H., & Zhang, Y. (2007) Geochimica et Cosmochimica Acta).


*** To find a decompression rate for a given set of input parameters user can vary conduit radius to achieve the desired final bubble number density. ***




Installation Guide:
This code runs on MATLAB 2020a. User need to have a valid MATLAB license and follow MATLAB installation guide. 

License: 
This code is under MIT license. 
https://opensource.org/licenses/MIT

Example:
Sample_Run.m provides an example of running the code.