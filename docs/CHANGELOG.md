# Changelog for crema  

## [Unreleased]

## [0.0.6] - 2023-06-02
### Added
- Added three peptide-level FDR procedures: psm-only, peptide-only, and
  psm-peptide
- Added two protein-level FDR procedures: max and sum
- Added mix-max
- Added support for MSAmanda, MSFragger, MSGF+, and Comet
- Added support for mzTab and pepXML file format

### Fixed
- Fixed a bug where newer versions of Pandas resulted in an error with peptide
  pairing.

## [0.0.5] - 2021-06-11  
### Added  
- This changelog.

### Changed  
- Refactored the structure of classes and functions for better organization and readability
- Updated and improved unit/system tests
- Huge speed increase in calculating confidence level estimates
