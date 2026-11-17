# Implementation notes

See the [README](./README.md) for setup and basic usage.

## Reader model and instrument detection

`BafReader` is the shared abstract interface and directory factory.
`BafReader.from_directory(...)` returns a `FticrBafReader` for FT-ICR data
or a `TofBafReader` for TOF data. Each concrete reader owns its calibration
cache and instrument-specific spectrum decoding: one FT-ICR calibration
versus TOF calibrations keyed by scan ID. Their `read_calibration()` methods
return `FticrCalibrationRecord` and `TofCalibrationRecord`, respectively.
Both record types are exported from `bruker_baf` and provide
`index_to_mz(index)`.

For direct construction from file paths, use the appropriate concrete
reader; `BafReader` itself is abstract. Explicit concrete readers reject
datasets from the other instrument family.

Instrument detection uses the null-terminated marker in the BAF preamble.
The model name comes from the acquisition header's string table.
`reader.instrument_info` exposes `instrument_type`, `model_name`,
`serial_number`, and `part_number`. Unknown or ambiguous instrument
markers raise `UnsupportedInstrumentError`.

## Acquisition header and scan selection

`reader.acquisition_header` includes the UTC start time and recorded offset,
analysis name, string table, variable definitions keyed by ID, and reference
calibration payloads. String-table entries are mapped to named optional
fields in `reader.acquisition_header.strings`:

- `source_path`, `instrument_name`, and `operator`
- `acquisition_software`, `software_version`, and `method_name`
- `dotted_version`, `role`, `sample_name`, and `sample_description`

`reader.scan_ids` lists the available scan IDs. Datasets with multiple scan
groups require an explicit `scan_id`; no scan is selected implicitly and
scans are not aggregated. `read_spectrum`, `read_calibration`, `read_maxima`,
`read_intensities`, and `read_metadata` accept scan selection.

FT-ICR defaults to profile mode; TOF defaults to centroid mode. FT-ICR
calibration still requires a single calibration record.

## FT-ICR calibration record payload

`FticrCalibrationRecord` represents FT-ICR (FTMS) mode-4 and mode-5
calibration. Records are stored as payloads inside `.baf_xtr` record
framing. The offsets below are relative to the calibration payload and
all numeric values are little-endian.

| Payload offset | Field | Type |
| --- | --- | --- |
| `0x00` | Declared payload size | `int32` |
| `0x04` | Header size | `int32` |
| `0x04 + header_size` | Block size | `int32` |
| Calculated offset | Calibration fields begin | See formula below |

The parser calculates the fields offset from the header:

```text
fields_offset = (4 + header_size) + (block_size - header_size)
```

For the observed layouts, `header_size` is 12. A block size of 60 gives a
fields offset of `0x40`; a block size of 64 gives `0x44`. Supported payload
sizes observed so far are 260, 300, and 308 bytes. The first word is the
declared payload size, not a format marker.

| Relative to fields offset | Field | Type |
| --- | --- | --- |
| `+8` | Mode | `int32` |
| `+12` | Alpha (mode 5) | `float64` |
| `+20` | Beta | `float64` |
| `+28` | High frequency | `float64` |
| `+40` | Size | `int32` |
| `+48` | Low frequency | `float64` |
| `+56` | Frequency width | `float64` |

## Native TOF spectrum support

### SQLite cache and calibration mapping

TOF reading uses an existing `analysis.sqlite` Baf2Sql cache alongside
`analysis.baf`, `analysis.baf_idx`, and `analysis.baf_xtr`. No Bruker DLL,
pyBaf2Sql, or additional Python dependency is needed. The cache is opened
read-only, never created or regenerated. Its recorded file sizes and
array identifiers are checked against the BAF files. Centroid arrays use
identity kinds `0x10`, `0x16`, and `0x11` for indices, m/z, and intensities.
Profile m/z and intensity identities use `0x43` and `0x42`. The lower 56
bits identify the corresponding BAF block offset.

Calibration lookup validates the identities for the indexed spectrum
blocks, without requiring centroid data. Profile-only scans use their
profile identities; centroid-only scans use their centroid identities.

`Spectra.TransformatorId` selects the scan's calibration from
`Transformators.Blob`. These blobs are also present in `.baf_xtr`; SQLite
provides the scan-to-calibration mapping, not an exclusive calibration
source. Native parsing currently supports the observed 112-byte mode-10
layout with zero higher-order correction. Reference transformators
(including mode 12) are not substituted for scan calibrations.

### Centroid layout and calibration equation

Centroid coordinates are float64 indices; intensities are float32 arrays
with their own payload offsets (not the FT-ICR float64 layout). For mode 10:

```text
time_ns = start_time + time_step * index - time_offset
time_ns = (1e6 / sqrt(quadratic_coefficient)) * sqrt(m)
          + mass_time_correction * m
m/z = m - mass_offset
```

The parser uses a numerically stable quadratic root.

### Profile layout and compression

The supported profile block has tag `0xBFA01001`, a 44-byte header, and
encoding field `0x00060102` at header offset 20. The word at offset 24 is
the payload byte count, which must equal the block size minus header size.
Other encodings are rejected, not interpreted as FT-ICR float arrays.

The observed compressed payload begins with little-endian uint32 words:
marker `0xEE77` and the decoded sample count. A fixed 30-byte codec
configuration follows, including its 26-byte settings section. Only the
verified configuration is accepted.

Data bits are read most-significant-bit first. Tokens use a unary prefix
of zero bits terminated by one:

| Leading zeros | Meaning |
| --- | --- |
| 0 | Repeat the previous intensity once |
| 1-7 | Magnitude bits followed by a sign bit |
| 8 | Full signed 32-bit delta |
| 9 | Unsigned 32-bit run length; repeat the previous intensity |

For prefixes 1-7 the magnitude bit widths are `4, 5, 6, 7, 9, 10, 11`;
their additive magnitude bases are `1, 17, 49, 113, 241, 753, 1777`.
A sign bit of one indicates a negative delta. The previous intensity
starts at zero and delta accumulation wraps as uint32, matching the
stored representation. The public API returns these values as floats,
like the other spectrum readers.

Runs must be nonempty and cannot exceed the remaining declared sample
count. After exactly that count, the decoder requires a 16-bit end marker
of `1`, zero padding to a 32-bit payload boundary, and a three-byte zero
block trailer. Truncation, invalid prefixes, overrunning runs, extra bytes,
and malformed markers or padding raise explicit errors.

Profile coordinates are consecutive sample indices starting at zero.
Each is converted to m/z using the selected scan's `TofCalibrationRecord`.
No centroid interpolation, binning, or scan aggregation is performed.

### Validation

Native output was checked against Bruker's Baf2Sql output for 209,735
centroids across two example datasets (97 scans), with maximum absolute
m/z error below `7e-13` Da and exact intensity matches. Native reading was
also exercised across all 17 cached example datasets (560 scans).

Profile intensities and m/z arrays were also compared in full with
Bruker's acquisition-calibration output on copied datasets. The first
three references cover 128 scans and 69,254,819 profile points: intensities
match exactly and maximum absolute m/z error is below `3e-12` Da.
All 560 cached example scans use the supported compression configuration.
The first and last profiles from each of the 17 datasets were decoded,
covering both observed profile lengths (618,833 and 297,678 samples).

### Limitations

Other TOF profile compression configurations, nonzero higher-order
corrections, and other calibration/centroid encodings raise
`NotImplementedError`.
The reader uses acquisition calibration from the cache, not later
recalibration in `calibration.sqlite` or `Calibrator.ami`, and does not
extract ion mobility.

Some TOF metadata value layouts (including block type 6 in the example
data) are also unsupported; a successful centroid read does not imply
that all metadata can be decoded.
