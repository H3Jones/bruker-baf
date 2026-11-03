# Bruker BAF reader

A Python reader for Bruker `.baf` mass-spectrometry data, derived from the
work of [rtms](https://nathanieltwarog.r-universe.dev/rtms).

## Supported data

- FT-ICR (FTMS) profile and centroid spectra.
- TIMS-TOF centroid spectra, using an existing `analysis.sqlite` Baf2Sql cache.
- Instrument identification and acquisition-header information.
- Explicit scan selection for datasets containing multiple scans.

Reading is native Python: no Bruker DLL or pyBaf2Sql dependency is required.
TOF profile spectra and ion-mobility extraction are not supported, and
some TOF metadata layouts cannot yet be decoded. TOF reading uses acquisition
calibration, not later recalibration.

## Setup

Requires Python 3.13 or newer. From a checkout, install the project and its
development dependencies with [uv](https://docs.astral.sh/uv/):

```powershell
uv sync
```

## Basic usage

Open a `.d` directory containing `.baf`, `.baf_idx`, and `.baf_xtr` files.
TOF datasets also require an existing `analysis.sqlite` cache; the reader
does not create or regenerate it.

```python
from bruker_baf import BafReader

reader = BafReader.from_directory(r"path\to\dataset.d")
print(reader.instrument_info)
print(reader.acquisition_header.started_utc)
print(reader.scan_ids)

scan_id = reader.scan_ids[0]
spectrum = reader.read_spectrum(scan_id=scan_id)
mz = spectrum.mass_to_charge
intensities = spectrum.intensities
```

The factory selects the appropriate instrument-specific reader. FT-ICR
defaults to profile mode; TOF defaults to centroid mode. Pass
`mode="centroid"` to request centroid spectra explicitly.

Scans are read individually, **not aggregated**. An explicit `scan_id` is
required when a dataset contains multiple scans; it may be omitted for a
single-scan dataset.

### Plot a spectrum

```powershell
uv run .\scripts\plot_spectrum.py "path\to\dataset.d" --plot
```

For a selected TOF scan:

```powershell
uv run .\scripts\plot_spectrum.py "path\to\tof.d" --scan-id 1 --mode centroid --plot
```

### Inspect acquisition strings

```powershell
uv run .\scripts\dump_acquisition_strings.py "path\to\dataset.d"
```

## Tests

```powershell
uv run pytest
```

## Implementation details

See [implementation notes](./implementation-notes.md) for reader design,
binary layouts, calibration equations, validation results, and detailed
limitations.
