# Surface variables verification

Verification of seasonal forecasting of surface variables (2m temperature, total precipitation and mean sea level pressure) for the Iberian Peninsula.

The calculation of the verification scores is based on the Copernicus Climate Change Service (C3S) [Seasonal Forecast Verification Tutorial](https://ecmwf-projects.github.io/copernicus-training-c3s/sf-verification.html).

## Index

* 1-Download_data.py: Python script for downloading seasonal forecasts and ERA5 data from the C3S Climate Data Store.
* 2-Postptocess.py: Python script for calclating verification scores.
* 3-Verification_plots.py: Python script for representing verification score maps for each forecasting system.
* 4-Multi-System_Verification_plots.py: Python script for representing verification score maps for all forecasting systems.
* 5-Plot_score-card.py: Python script for representing verification score-cards, with all systems and variables.
* 0-driver.sh: Bash script for automatically running all scripts, for all systems and variables.
* modesl.csv: List of all most recent forecasting systems available at the C3S Climate Data Store.
