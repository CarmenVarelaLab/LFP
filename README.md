# LFP

Group of matlab Functions for LFP analysis.

General guidelines:  

Do not use analysis code as a black box. Look at the raw data and processed data. 
Use scripts as a starting point. E.g., event detection scripts may require parameters to be adjusted for different brain regions, subjects, etc.

Detect Sleep:  Estimates putative quiet wakefulness/sleep states given:
    delta and ripple filtered LFPs (from Cortex and hippocampus respectively) if both available, or from Cx-delta alone
