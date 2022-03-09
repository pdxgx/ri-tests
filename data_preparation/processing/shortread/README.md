# Description

Scripts for running five short-read retained intron-detection (SR RI-detection) tools, formatting results as GenomicRanges, calculating the length-weighted medians (LWMs) of SR tool outputs, identifying called RIs from tool-specific RI confidence filters, and harmonizing short- and long-read expression together.

# Script descriptions
* `run_superintronic` -- Runs the SR RI-detection tool `superintronic`.
* `run_kma` -- Runs the SR RI-detection tool Keep Me Around (`KMA`).
* `run_irfinders` -- Runs the SR RI-detection tool `IRFinder-S`.
* `run_iread_ips` -- Runs the SR RI-detection tool `iREAD`.
* `run_interest` -- Runs the SR RI-detection tool `IntEREst`.
* `mpileup_server` -- Calculates median gene coverage using pileup strategy with `samtools`.
* `iread_make_results_table` -- Formats results output from `iREAD` for further analyses.
* `harmonized_lwm` -- Calculates harmonized length-weighted medians (LWMs) for five SR IR-detection tools using long-read intron ranges.
* `format_longread` -- Formats long-read results output for harmonization with short-read results.
* `format_superintronic_results` -- Formats results from `superintronic` as GenomicRanges for harmonization/LWM calculations.
* `format_kma_results` -- Formats results from `KMA` for as GenomicRanges for harmonization/LWM calculations.
* `format_irfinders` -- Formats results from `IRFinder-S` as GenomicRanges for harmonization/LWM calculations. 
* `format_iread_results` -- Formats results from `iREAD` as GenomicRanges for harmonization/LWM calculations.
* `format_interest_results` -- Formats results from `IntEREst` as GenomicRanges for harmonization/LWM calculations.