# Changes in locomotion and spreading depolarizations after photothrombosis
Analysis scripts related to the manuscript "Unilateral hippocampal stroke in freely behaving mice reveals sex differences in contralesional spreading depolarization and related behaviour"

## Prerequisites
- [MATLAB][MATLAB] (last tested with R2023a)
- [FPA library][FPA]

## Installation
- Install [MATLAB][MATLAB] with the following toolboxes:
	- Curve Fitting Toolbox
	- Signal Processing Toolbox
	- Statistics and Machine Learning Toolbox
- Download and extract the [FPA library][FPA] to the `Documents/MATLAB` folder
- Download and extract these scripts to the `Documents/MATLAB` folder.

## Usage
- Go to the FPA folder and run `startup.m` to add FPA dependencies to the MATLAB search path.
- Edit `analysis1.m` to match your experimental configuration (e.g. path to files and column numbers, etc) then run to pre-process fibre-photometry data and get behavior triggered data.
- Run `analysis2.m` to detect peak derivatives in ipsilesional GCaMP6f related to stroke.

## Changelog
See [Changelog](CHANGELOG.md)

[FPA]: https://github.com/leomol/FPA
[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md