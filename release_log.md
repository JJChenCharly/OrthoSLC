# Release Log OrthoSLC (0.1.1)

`0.1Alpha` -> `0.1Beta`:<br>
* allow users to set `bin level`, in steps using hash binning, with a linear style instead of exponatial style.<br>

* **Fixed**: fixed output path concatenation failure in Step 2.

`0.1Beta` -> `0.1`:<br>
* applied pre-clustering before BLAST which significantly reduce time usage especially when phylogenetically close genomes paticipated analysis.
    * 500 <i>Listeria monocytogenes</i> in `0.1Beta` -> ~7 hours, `0.1` -> ~22 mins
    * 1150 <i>E. coli</i> `0.1Beta` -> ~2.7 days, `0.1` -> less than 4 hours
    
* In `Step5_reprocal_blast`, we provide optinal memory efficient mode. It is a simple set up of multithreading when running BLAST. Since after pre clustering, some non-redundant genomes may become very small in size, which make them less efficient to BLAST with multiple threads. With memory efficient mode `off`, program will assign only 1 thread on each BLAST task so that there will be no thread waiting for database formatting and pre-processing. It can significantly increase the BLAST efficiency but with some more consumed memory.

* we provide **no lock mode** in all steps that apply hash binning to speed up the process. We allow users to turn off mutex lock which is to safely write into files when multi-threading. In ours tests, program can generate files without data corruption when multi-threading with no lock (data corruption were rarely observed, the possiblity of data corruption may vary between computation platform).

* **Fixed**: fixed the bug happening when writing final clusteres of FASTA files.

`0.1` -> current version `0.1.1`:<br>
* smarter thread assignment in `Step5_reprocal_blast`.
* more complete file and directory path check in callable binary.
