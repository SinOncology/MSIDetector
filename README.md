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

Output example
------------
Below is the output of MSIDetector:

<pre>
loci	sampleModelValue	sampleModelPvalue	modelRepeatCountIndexes	modelMean	modelSd	modelValueThreshold	modelPvalueThreshold	sampleLociLabel
chr11-102193508-102193534-NR-27	0.990885	0.0381868	18,19,20,21,22,23,24,25,26,27,28	0.996227	0.0024656	0	3.21166e-17	stable
chr11-108114661-108114676-T15-1	0.982578	1.01231e-12	13,14,15,16,17	0.997162	0.00199567	0	6.64371e-28	stable
chr11-108141955-108141970-T15-2	0.883752	0.398896	14,15,16	0.883432	0.0208629	0	0.00443154	stable
chr11-108188266-108188279-T13	0.997283	0.321196	11,12,13,14,15	0.995515	0.00268452	0	2.22712e-08	stable
chr11-125490765-125490786-NR-22	0.98727	0.0886214	17,18,19,20,21,22,23,24	0.994232	0.00401337	0	1.68345e-06	stable
chr14-23652346-23652367-NR-21	0.993256	0.391769	18,19,20,21,22,23,24,25,26,27	0.987344	0.0310385	0	0.36861	stable
chr16-56718015-56718035-MT1X	0.992055	0.287307	15,16,17,18,19,20,21,22,23	0.994229	0.00268405	0	2.85343e-12	stable
chr17-29508819-29508835-T16	0.980926	0.137362	14,15,16,17,18,19	0.988553	0.00522284	0	5.87209e-13	stable
chr2-39536689-39536716-MONO-27	0.995062	0.298611	20,21,22,23,24,25,26,27,28,29,30,31,32,33	0.997146	0.00273892	0	2.53284e-15	stable
chr2-47641559-47641586-BAT26	0.997773	0.3186	16,17,18,19,20,21,22,23,24,25,26,27,28,29,30	0.998844	0.00159695	0	1.91513e-17	stable
chr2-95849361-95849384-NR24	0.996217	0.397084	17,18,19,20,21,22,23,24,25,26	0.995972	0.00253186	0	9.10046e-10	stable
chr4-153268227-153268241-A14	0.983884	0.293527	12,13,14,15,16	0.986729	0.00363149	0	1.4859e-06	stable
chr4-55598211-55598236-BAT25	0.985556	0.371132	20,21,22,23,24,25,26,27,28	0.989925	0.0114944	0	0.000957318	stable
chr6-117725383-117725398-A15	0.941828	0.106558	14,15,16,17	0.956682	0.00914136	0	0.000983155	stable
chr7-116409675-116409690-T15-3	0.0185615	0.0520716	9,10,11	0.00872715	0.00487327	0	1.63607e-06	stable
chr7-143003342-143003367-CAT25	0.964286	0.000236134	20,21,22,23,24,25,26,27,28	0.987817	0.00610337	0	1.48672e-06	stable
chr9-135773000-135773018-A18	0.993711	0.243921	15,16,17,18,19,20,21,22	0.988359	0.00539492	0	0.00121948	stable
<br>
#Judge Criteria:
#  Unstable Rate is in [0,0.2): MSI-Stable
#  Unstable Rate is in [0.2,0.4): MSI-Low
#  Unstable Rate is in [0.4,1]: MSI-High
#Conclusion:
#  Unstable loci count is: 0, Total Loci Count is: 17
#  The unstable rate is: 0
#  The sample MSI status is: MSI-Stable
</pre>

Citation
--------

Jin et al.
[Detection of Microsatellite Instability Status from Targeted Sequencing Data.]
Bioinformatics, 2019 submitted.


License
-------
MSIDetector is distributed under the GPLv3. Consult the accompanying [LICENSE](https://github.com/SinOncology/MSIDetector/LICENSE) file for more details.
