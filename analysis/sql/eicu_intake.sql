SELECT patientunitstayid
, sum(cellvaluenumeric) AS sum_intake
FROM `physionet-data.eicu_crd.intakeoutput` 
WHERE LOWER(cellpath) LIKE '%intake (ml)%' 
AND intakeOutputOffset BETWEEN -6*60 AND 6*60
GROUP BY patientunitstayid
ORDER BY patientunitstayid
LIMIT 1000
