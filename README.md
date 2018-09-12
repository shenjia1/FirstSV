# FirstSV

This tool is used to detect structural variants from a BAM file using local assembly.

# All options

```shell
usage: FirstSV V 0.02 [-h] --bamfile bamfile --outputdir outdir --sampleid
                      sampleid --configinfo config_file
                      [--targetregion targetregion]
                      [--bpthreshold bpthreshold] [--spthreshold spthreshold]
                      [--mpthreshold mpthreshold] [--hard_filter hard_filter]

This tool is used to detect structural variants from a BAM file using local
assembly.

optional arguments:
  -h, --help            show this help message and exit
  --bamfile bamfile, -b bamfile
                        BAM file created by BWA or other alignment software.
  --outputdir outdir, -o outdir
                        Output direction for results.
  --sampleid sampleid, -i sampleid
                        Sample ID.
  --configinfo config_file, -c config_file
                        Config file direction.
  --targetregion targetregion, -r targetregion
                        Detect Region.
  --bpthreshold bpthreshold, -p bpthreshold
                        Breakpoint Mini Distance.
  --spthreshold spthreshold, -s spthreshold
                        Breakpoint Mini Support Counts.
  --mpthreshold mpthreshold, -m mpthreshold
                        Mini Mapping Quality.
  --hard_filter hard_filter, -f hard_filter
                        Weather to using hard filter.(EXPERIMENTAL)

```

# Simple usage
```shell
python FirstSV.py -c config.txt -o /PATH_OF_OUTPUT -i SAMPLE_ID -b BAM_FILE
```

# Dependencies
* Standalone BLAT v. 36x1
* samtools-1.3
* bwa-0.7.13

#  Acknowledge
fermi_las was based on the [https://github.com/lh3/fermikit](fermikit) created by LiHeng.

