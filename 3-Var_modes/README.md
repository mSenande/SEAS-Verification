# Climate variability modes verification

Verification of seasonal forecasting of climate variability modes (North Atlantic Oscillation, East Atlantic, East Atlantic / Western Russia and Scandinavian Pattern).

The calculation of the verification scores is based on the Copernicus Climate Change Service (C3S) [Seasonal Forecast Verification Tutorial](https://ecmwf-projects.github.io/copernicus-training-c3s/sf-verification.html).

## Index

* 1-Box_calc: Verification of seasonal forecasting of North Atlantic Oscillation calculated as the difference between the normalized anomaly of sea level pressure between a South-Atlantic and a North-Atlantic lat-lon box. The methodology is based on Baker et al., 2018.
* 2-EOFs_calc: Verification of seasonal forecasting of four main variability modes of each season calculated with an Emprirical Orthogonal Function analysis. The methodology is based on Lledó et al., 2020.
* 2b-EOFs_calc: Verification of seasonal forecasting of four specific variability modes (North Atlantic Oscillation, East Atlantic, East Atlantic / Western Russia and Scandinavian Pattern) calculated with an Emprirical Orthogonal Function analysis of the DJF season.
* 2c-EOFs_calc: Verification of seasonal forecasting the four main variability modes, calculated with an Emprirical Orthogonal Function analysis (using the mean sea level pressure instead of the 500hPa geopotential height).
* 3-Common_plots: Scripts for representing verification score-cards for a particular season, with all variability modes and both of the above-mentioned methodologies.

## References

Baker, L. H., Shaffrey, L. C., Sutton, R. T., Weisheimer, A., & Scaife, A. A. (2018). An intercomparison of skill and overconfidence/underconfidence of the wintertime North Atlantic Oscillation in multimodel seasonal forecasts. *Geophysical Research Letters*, 45(15), 7808-7817.

Lledó, L., Cionni, I., Torralba, V., Bretonniere, P. A., & Samso, M. (2020). Seasonal prediction of Euro-Atlantic teleconnections from multiple systems. *Environmental Research Letters*, 15(7), 074009.
