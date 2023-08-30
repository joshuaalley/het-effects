* replicate all analyses
version 16.1

* Rename and Unzip replication archive
if "`c(os)'" == "Windows" {
   shell ren TomzWeeks-Alliances-ISQ-Files TomzWeeks-Alliances-ISQ-Files.zip
}
else {
   shell mv TomzWeeks-Alliances-ISQ-Files TomzWeeks-Alliances-ISQ-Files.zip
}
unzipfile TomzWeeks-Alliances-ISQ-Files.zip

* Study 1: 2017-04-YouGov
cd 2017-04-YouGov/Code
do 2017-04-YouGov-runner.do

* Study 2: 2017-12-Lucid
cd ../../2017-12-Lucid/Code
do 2017-12-Lucid-runner.do

* Study 3: 2020-09-Lucid
cd ../../2020-09-Lucid/Code
do 2020-09-Lucid-knowledge.do
