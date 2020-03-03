# Striatal modelling
Striatal microcircuit generation and data analysis

## Overview
This code generates a 'physical' striatal microcircuit based on the parameters described in [Humphries *et al.*, 2010](https://doi.org/10.1371/journal.pcbi.1001011) and outputs neural connection lists in CSV format, suitable for importing into the [SpineCreator](http://spineml.github.io/spinecreator/) modelling tool.

## Requirements
Several additional MATLAB functions are required:
* [cell2csv](https://mathworks.com/matlabcentral/fileexchange/7601-cell2csv)
* [Chronux](http://chronux.org) for analysing spiking oscillations
* [geom3d](https://mathworks.com/matlabcentral/fileexchange/24484-geom3d)
* [ls2](https://mathworks.com/matlabcentral/fileexchange/40042-recursive-directory-search)
* [struct2csv](https://mathworks.com/matlabcentral/fileexchange/34889-struct2csv)
* [xml2struct](https://mathworks.com/matlabcentral/fileexchange/28518-xml2struct)
