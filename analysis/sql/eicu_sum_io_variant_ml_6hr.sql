SELECT patientunitstayid
, sum(CASE when LOWER(cellpath) LIKE '%intake (ml)%' 
  AND intakeOutputOffset BETWEEN -6*60 AND 6*60
  Then cellvaluenumeric Else 0 End) AS intake_ml_6hr
, sum(CASE when LOWER(cellpath) LIKE '%output (ml)%' 
AND intakeOutputOffset BETWEEN -6*60 AND 6*60
  Then cellvaluenumeric Else 0 End) AS output_ml_6hr
, sum(CASE when LOWER(cellpath) LIKE '%intake (ml)%' 
  AND intakeOutputOffset BETWEEN -6*60 AND 6*60
  Then cellvaluenumeric Else 0 End) 
- sum(CASE when LOWER(cellpath) LIKE '%output (ml)%' 
  AND intakeOutputOffset BETWEEN -6*60 AND 6*60
  Then cellvaluenumeric Else 0 End) AS io_variant_ml_6hr
FROM `physionet-data.eicu_crd.intakeoutput` 
GROUP BY patientunitstayid
ORDER BY patientunitstayid
LIMIT 1000
