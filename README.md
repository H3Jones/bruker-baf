# Bruker .baf file parser

Rough python implementation derived from the work of [rtms](https://nathanieltwarog.r-universe.dev/rtms)


## Calibration record payload

Calibration records are stored as payloads inside the `.baf_xtr` record
framing. The offsets below are relative to the calibration payload and all
numeric values are little-endian.

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