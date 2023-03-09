This archive contain the input and the output of the instances solved in Section 7 of the paper The service points' location and capacity problem, Raviv (2023).

The demand points of the three cities are stored in the files <city_name>-dp.csv.  For each point we store the lat, lon an poulation size.

The SP candidate locations are stored in files <city_name>-sp-map-<model (cover/PWL)>-<demand (high/low)>.csv
Low demand refer to expcted demand of 2% of the population of each demand point and high demand to 4%.
The files contains the following columns
LAT,LON - of the candidate SP locations
Type - the type of the SP opened in the location (0-no, 1 - SP with capcity 30,  2-60, 3-100, 4-150)


