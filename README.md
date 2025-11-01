# ONTstats
## Summary
This repository contains scripts to calculate statistics of the nanopore sequencing data produced in the [PuckerLab](https://www.izmb.uni-bonn.de/en/pbb) on a P2 Solo.

## Usage
Python3 and the library matplotlib are required to run this analysis.

```
  python3 evaluate_ONT_P2_runs.py --in <INPUT_FILE> --out <OUTPUT_FOLDER>
  --in           STR     Input file name
  --out          STR     Output directory
  --nomasking            Inactivates the masking of flow cell IDs
```

Details about the arguments:
`--in` specifies the input file which contains the details about each nanopore sequencing run.

`--out` specifies the output folder which will contain the statistic figures. If the folder does not exist, it will be created.

`--nomasking` deactivates the masking of flow cell IDs. Default: replacement of IDs by running numbers.

## Example data
The example data set allows to test the Python script. Flow cell IDs and values have been altered in the exmaple data set.

## References

