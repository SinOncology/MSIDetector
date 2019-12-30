`<p align="center">
  <a href="http://www.sinotechgenomics.com">
    <img height="70" src="http://www.sinotechgenomics.com/Upload/0/WebsiteLogo/WebsiteLogo_20170620092534731.png">
  </a>
  <h1 align="center">BreakID</h1>
</p>


MSIDetector: Detection of Microsatellite Instability Status from Targeted Sequencing Data without Matched Normal
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

`bin/MSIDetector -i input_bam -o out_prefix  -l /path/to/MSI_Loci_File -ref /path/to/reference_baseline_file`


Citation
--------

Jin et al.
[Detection of Microsatellite Instability Status from Targeted Sequencing Data.]
Bioinformatics, 2019 submitted.


License
-------
MSIDetector is distributed under the GPLv3. Consult the accompanying [LICENSE](https://github.com/SinOncology/MSIDetector/LICENSE) file for more details.
