# Changelog

Human-readable summary of major changes and features/data that are added or removed.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
### Added
### Removed


## [0.2.0] - 2022-02-09
Reviewed for issues with the code, removed dangling branches, added python help menus and dependencies in setup.py.
### Changed
- Formalized regression tests, added annotations to plots and explicit checks.  Restructed these directories.
### Added
- Added `install_requires` to `setup,py`, changed install instructions.
- Added argparse to `do_map`, `make_ref`, and `smooth_ref`, as well as help functions.  Note that this changes the API.
### Removed
- Various shells like `run_map.sh`, etc, which remove some of the indirection.

