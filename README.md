`<p align="center">

  <h1 align="center">MSIDetector</h1>
</p>


MSIDetector: Detection of Microsatellite Instability Status from Targeted Sequencing Data


Prerequisites for MSIDetector
----------------
cmake version above 2.8.2

Installing MSIDetector
----------------

`git clone https://github.com/SinOncology/MSIDetector.git`

`cd MSIDetector/`

`sh install.sh`

Running MSIDetector
--------------------------
MSIDetector needs a sorted, indexed and duplicate marked tumor bam file, a MSI loci file and pre-build reference baseline file. The output is in plain txt format.

`bin/MSIDetector -i input_bam -o out_prefix  -l /path/to/MS_Loci_File -ref /path/to/reference_baseline_file`

Example MS loci file :
https://github.com/SinOncology/MSIDetector/blob/master/baseline_data/msi_loci.bed

Example reference baseline file
https://github.com/SinOncology/MSIDetector/blob/master/baseline_data/msi_loci_baseline.txt


Citation
--------

Jin et al.
[Detection of Microsatellite Instability Status from Targeted Sequencing Data.]
Bioinformatics, 2019 submitted.


License
-------
MSIDetector is distributed under the GPLv3. Consult the accompanying [LICENSE](https://github.com/SinOncology/MSIDetector/LICENSE) file for more details.
