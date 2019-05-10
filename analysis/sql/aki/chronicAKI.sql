-- Chronic patients with AKI receiving rrt prior to ICU admission
SELECT
  DISTINCT patientunitstayid
FROM
  `physionet-data.eicu_crd.treatment`
WHERE
  LOWER(treatmentstring) LIKE '%rrt%'
  OR LOWER(treatmentstring) LIKE '%dialysis%'
  OR LOWER(treatmentstring) LIKE '%ultrafiltration%'
  OR LOWER(treatmentstring) LIKE '%cavhd%' 
  OR LOWER(treatmentstring) LIKE '%cvvh%' 
  OR LOWER(treatmentstring) LIKE '%sled%'
  AND 
  LOWER(treatmentstring) LIKE '%chronic%'
UNION  DISTINCT

  SELECT
  DISTINCT patientunitstayid
FROM 
  `physionet-data.eicu_crd.apacheapsvar`
WHERE
  dialysis = 1 -- chronic dialysis prior to hospital adm