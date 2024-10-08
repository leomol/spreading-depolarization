## Description
Analysis scripts related to the manuscript "Changes in open field activity, bilateral local field potential (ie. epileptiform activity) and bilateral GCaMP photometry (ie. spreading depolarizations) after unilateral photothrombosis"

## Prerequisites
- [MATLAB][MATLAB] (last tested with R2023b)
- [FPA library][FPA]

## Installation
- Install [MATLAB][MATLAB]
- Download and extract the [FPA library][FPA] to the `Documents/MATLAB` folder (see [FPA's installation instructions][FPA] for the required toolboxes)
- Download and extract the [spreading-depolarization][project] scripts to the `Documents/MATLAB` folder.

## Usage
- Go to the `FPA` folder and run `startup.m` to add FPA dependencies to the MATLAB search path.
- Go to the `spreading-depolarization` folder and run `startup.m` to download example data. Remove `startup.m` and/or the data if you do not require these data.
- Edit `analysis1.m` to match your experimental configuration (e.g. path to files and column numbers, etc) then run to pre-process fibre-photometry data to get behavior triggered data.
- Run `analysis2.m` to detect peak derivatives in ipsilesional GCaMP6f related to stroke.
- Run `analysis3.m` to detect epileptiform activity in bilateral local field potential recordings relative to simultaneous bilateral fibre photometry recordings.

## Changelog
See [Changelog](CHANGELOG.md)

[FPA]: https://github.com/leomol/FPA
[project]: https://github.com/leomol/spreading-depolarization
[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md