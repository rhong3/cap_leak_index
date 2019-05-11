-- ------------------------------------------------------------------
-- Title: Patients intake (ml)
-- Notes: cap_leak_index/analysis/sql/eicu_intake.sql 
--        cap_leak_index, 20190511 NYU Datathon
--        eICU Collaborative Research Database v2.0.
-- ------------------------------------------------------------------
SELECT patientunitstayid
, sum(cellvaluenumeric) AS sum_intake
FROM `physionet-data.eicu_crd.intakeoutput` 
WHERE LOWER(cellpath) LIKE '%intake (ml)%' 
AND intakeOutputOffset BETWEEN -6*60 AND 6*60
GROUP BY patientunitstayid
ORDER BY patientunitstayid
LIMIT 1000
