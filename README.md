# spectral 1D extraction

## Requirements:
* CarPy

## Usage examples:

Run the 1D extraction
`extract_1dspec.py -smf 021953.SMF -dir 021953`

Run the spectral viewer
`g_spec.py -smf 021953.SMF -dir 021953`
`g_spec.py -smf 021953.SMF -dir 021953`--trace 021953.trace

## Keyboard commands for g_spec.py:

| key | description |
| --- | ---|
| nN | next/previous |
| hjkl | left/down/up/right |
| sdbf | sharper/duller/brighter/fainter |
| g | go to frame? |
| r | reset display |
| z | turn on/off redshift lines display |
| Z | enter redshift? |
| H | plot redshift histogram |
| B | boxcar width? |
| 012345 | redshift quality |
| [] | increment/decrement redshift by 0.01 |
| {} | increment/decrement redshift by 0.001 |
| C | enter comment? |
| p | print all object information |
| w | enter wavelength of line cursor is over? |
| T | display trace (if loaded) |
| Q | save and quit |




