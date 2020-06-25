# spectral 1D extraction

## Requirements:
* CarPy 

## Usage examples:

Run the 1D extraction

`extract_1dspec.py -smf 021953.SMF -dir 021953`

Run the spectral viewer

`g_spec.py -smf 021953.SMF -dir 021953`

`g_spec.py -smf 021953.SMF -dir 021953 --trace 021953.trace`

## Keyboard commands for g_spec.py:

| key | description |
| --- | ---|
| nN | next/previous |
| hjkl | left/down/up/right |
| sdbf | sharper/duller/brighter/fainter |
| g | go to frame? |
| r | reset display |
| R | rebinning (x,y)? |
| z | turn on/off redshift lines display |
| Z | enter redshift? |
| H | plot redshift histogram |
| B | boxcar width? |
| [] | increment/decrement redshift by 0.1 |
| ;' | increment/decrement redshift by 0.01 |
| ,. | increment/decrement redshift by 0.001 |
| {} | increment/decrement redshift by 0.0001 |
| C | enter comment? |
| S | sum emission line flux between two points |
| p | print all object information |
| w | enter wavelength of line cursor is over? |
| T | display trace (if loaded) |
| q | z-quality? |
| X | create ds9 region file for redshifts |
| Q | save and quit |




