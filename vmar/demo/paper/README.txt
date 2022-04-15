Code to reproduce figures for Geophysics paper:

files:
direct.scr  ExtractingTimeDifferences.py  isolation.scr  model.scr  pickedwindows.npz  README.txt  remove_direct.scr  shots_slurm.scr  supython.py  Td_calc.scr

1. Marchenko based isolation:
   a. model.scr --> create baseline and monitor model
   b. shots_slurm.scr --> create baseline and monitor reflection response
   c. direct.scr --> create direct arrival for all shots
   d. remover_direct.scr --> removes the direct arrival from the reflection responses
   e. Td_calc.scr --> create initial focus function and two-way travel times for mute windows
   f. isolation.scr --> isolate the target response from the baseline/monitor reflection response
     jobtime ~1h.

2. Extracting time-differences:
   - ExtractingTimeDifferences.py --> python code to get the time-difference from the zero-offset sections (tested on Windows and Unix)
   - Close the first picking window that pops up to continue to code, the picking used in the paper is provided by: "pickedwindows.npz", but can interactively be changed in this first figure
   - note that this requires supython to read su files into python as Numpy arrays (see https://github.com/Jbrackenhoff/SUpython)

