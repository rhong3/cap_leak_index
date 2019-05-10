WITH
first_hct_6hrs AS (
  SELECT
  patientunitstayid,
  hematocrit AS first_hct_6hrs,
  ROW_NUMBER() OVER (PARTITION BY patientunitstayid ORDER BY chartoffset ASC) AS position
  FROM
  `physionet-data.eicu_crd_derived.pivoted_lab` pivoted_lab
  WHERE
  chartoffset BETWEEN -6*60
  AND 6*60 ),
mean_hct_24_36hrs AS (
  SELECT
  patientunitstayid,
  ROUND( AVG (CASE
              WHEN chartoffset BETWEEN 24*60 AND 36*60 AND hematocrit IS NOT NULL THEN hematocrit
              END
  ),2) AS mean_hct_24_36hrs
  FROM
  `physionet-data.eicu_crd_derived.pivoted_lab`
  GROUP BY
  patientunitstayid )
SELECT
first_hct_6hrs.patientunitstayid,
first_hct_6hrs.first_hct_6hrs,
mean_hct_24_36hrs.mean_hct_24_36hrs
FROM
first_hct_6hrs
INNER JOIN
mean_hct_24_36hrs
USING
(patientunitstayid)
WHERE
position = 1
and mean_hct_24_36hrs IS NOT NULL