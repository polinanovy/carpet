## Searching for event excess from Cygnus Cocoon region

**cygnus.py** - the main program
Calculate angle between Carpet's events and Cygnus Cocoon (RA = 307.18, DEC = 41.31)
Calculate Signal-to-Noise ratio

Photon-like events were selected in two ways: 
1. $n_{\mu}$ = 0
2. $(n_{\mu} + 0.1) / N_e > 10^{-5.90688}$

Hence, we have 3 types for each result file, for instance:
- **CC_match.txt** - list of all Carpet's events matched with Cygnus Cocoon
- **CC_match_photons_nmu0.txt** - list of photon-like events selected in the first way matched with Cygnus Cocoon
- **CC_match_photons.txt** - list of photon-like events selected in the second way matched with Cygnus Cocoon
