 # Notes:

## Domain

- Pab1 - domain is [(125, 200)]
- Bgl3 - no domain because PBD starts as same place as protein sequence (Boolean of PDB lines up with protein without change)
- Ube4b - domain is [(1071, 1173)], but position in Ube4b


## STRIDE

Documentation:

       H	    Alpha helix
	   G	    3-10 helix (variation of alpha helix)
	   I	    PI-helix
	   E	    Extended conformation (or strand - beta strands?)
	   B or	b   Isolated bridge (one hydrogen bond of the type found in a beta-sheet)
	   T	    Turn
	   C	    Coil (none of the above) (absence of regular secondary structure)


 - currently not including **C** or **T**


Meeting Notes 6/23 (next meeting July 5th)
- Limiting the number of mutations in sec. str. vs. out of sec. str. (Friday)
  - changing the blocks to be equal
- running again with the limiting number (do it again with the chunks?) (Monday)
- running with alpha helices versus beta sheets (Tuesday)
- finding ~10 proteins with functional data (Tuesday/Wednesday)
- run all 10 (Wednesday/Thursday)
- fraction of sec. str. vs. learnability


Questions
- why is pab1 consistently higher in terms of secondary structure?


## Choosing Other proteins
- 3 from curated database:
   - nuclease - (ddG 1282 vals)
   - lysozyme - (ddG 1376 vals)
   - protein G - (ddg 1221 vals)
