# Azimuth vs CellTypist vs PopV

Code:
azimuth_prediction_pipeline.R - Transfers labels onto test RDS-dataset based on reference dataset anchors. Output can be configured to be CSV/RDS file.
Azimuth_convert_Rds_to_h5ad.R - Script takes in command-line parameters to convert RDS into h5ad. Need this for future downstream analytics.

Datsets:
LungMap.h5ad- Source data (https://celltype.info/CAPinitialRelease/LungMAP-Human-data-from-a-broad-age-healthy-donor-group/3)
PopV_Reference_Lung.h5ad- Training dataset for PopV model (https://www.dropbox.com/s/2kuzdamjevev2ci/Lung.h5ad?dl=1)
