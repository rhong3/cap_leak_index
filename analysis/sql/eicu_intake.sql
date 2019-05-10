SELECT patientunitstayid
, sum(cellvaluenumeric) AS sum_intake
FROM `physionet-data.eicu_crd.intakeoutput` 
WHERE LOWER(cellpath) LIKE '%intake (ml)%' 
GROUP BY patientunitstayid
LIMIT 1000
