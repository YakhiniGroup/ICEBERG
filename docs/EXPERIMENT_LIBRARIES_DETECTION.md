# Step Title


Each paired-end reads, r1 and r2, from R1.fastq and R2.fastq respectively, is from one of the following libraries:<br>
Forward library: paired-end reads where part of the experiment tag reverse complement is at the start of r2.<br>
Reverse library: paired-end reads where part of the experiment tag is at the start of r2.<br>

This happens due the PCR rounds for amplify sites containing the experiment injected tag, The PCR primers are - <br>
Forward primer: tag_reverse_complement\[12:] (22 nt). <br>
Reverse primer: tag\[6:] (28 nt). <br>
<br><br>
 
 
The library of each paired-end reads is detected by aligning the libraries PCR-primers to 
the start of r2, r2\[:22] or r2\[:28] respectively, the library of the primer with the maximum alignment-score 
is the one detected for both r1 and r2.


The detected PCR primer is appended to the (identical) umi section in r1 and r2 fastq headers.


These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Prerequisites

The things you need before installing the software.

* You need this
* And you need this
* Oh, and don't forget this

## Usage

A few examples of useful commands and/or tasks.

```
First example
Second example
And keep this in mind
```
