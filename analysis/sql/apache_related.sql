-- ------------------------------------------------------------------
-- Title: Select patients from apachepatientresult 
-- Notes: cap_leak_index/analysis/sql/apache_related.sql
--        cap_leak_index, 20190511 NYU Datathon
--        eICU Collaborative Research Database v2.0.
-- ------------------------------------------------------------------
SELECT patientunitstayid, apachescore, actualicumortality, 
actualhospitalmortality, unabridgedunitlos, unabridgedhosplos, unabridgedactualventdays 
FROM `physionet-data.eicu_crd.apachepatientresult`
