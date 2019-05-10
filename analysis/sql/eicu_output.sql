SELECT patientunitstayid
, sum(cellvaluenumeric) AS sum_output
FROM `physionet-data.eicu_crd.intakeoutput` 
WHERE LOWER(cellpath) LIKE '%output (ml)%' 
GROUP BY patientunitstayid
LIMIT 1000
