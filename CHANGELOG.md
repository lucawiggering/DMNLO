# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2023-11-16

### Changed

- The mass mixing matrices are now diagonal for the first two generations. 
- The charm mass is set to zero in 'DMNLO_MSSMinit.F'.
- The W- and Z-widths are now hard coded according to their PDG 2022 values.
- Before m_t^{DR,MSSM} was computed using the Mteff function from MicrOmegas. Now it is computed with the mOS2mMS function from RunDec.
- Support for the photon final state in NeuQ2qx was removed because the cancellation of the IR divergences from the soft photon require the inclusion of O(alpha) virtual corrections.


