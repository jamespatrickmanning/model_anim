Originally Lei Zhao's code from 2019 to animate DOPPIO model temps

Modified by JiM in Feb 2020 to allow user specified "area" and hourly images by default. 

Modified by Mingchao in Feb 2020 select good temperature contour levels and added GOMOFS

Modified by JiM in Summer 2020 with Jack Polentes help

Modified by JiM in Fall 2020 by adding drifter tracks and, in late Oct 2020, NCEP wind vector

Modified by JiM in Feb 2021 to have user select between DOPPIO,  GOMOFS, FVCOM w/"make_model_gif.py"

Stored in Github as "model_anim" repository

See hardcodes at top of code where there is, for example, "start_date" and "ndays" 
Note: You may need to adjust the "clevs" to get good range of colors.
Note: You may want to adjust the Basemap resolution to "c" for crude in experiments.
Note: You may need to "conda install -c conda-forge basemap-data-hires" in order to use the higher resolution coastlines

Requires the module "zlconversions.py" which have many functions.
