WITH
  blood AS(
  SELECT
    patientunitstayid
  FROM
    `physionet-data.eicu_crd.intakeoutput`
  WHERE
    intakeoutputentryoffset BETWEEN -6*60
    AND 36*60
    AND LOWER(cellpath) LIKE '%blood%' ),
  others AS (
  SELECT
    patientunitstayid,
    SUM(CASE
        WHEN LOWER(cellpath) LIKE '%output (ml)%' AND LOWER(cellpath) LIKE '%gastric%' AND intakeOutputOffset BETWEEN -6*60 AND 6*60 THEN cellvaluenumeric
      ELSE
      0
    END
      ) AS gastric,
    SUM(CASE
        WHEN LOWER(cellpath) LIKE '%output (ml)%' AND LOWER(cellpath) LIKE '%stool%' AND intakeOutputOffset BETWEEN -6*60 AND 6*60 THEN cellvaluenumeric
      ELSE
      0
    END
      ) AS stool,
    SUM(CASE
        WHEN LOWER(cellpath) LIKE '%output (ml)%' AND LOWER(cellpath) LIKE '%emesis%' AND intakeOutputOffset BETWEEN -6*60 AND 6*60 THEN cellvaluenumeric
      ELSE
      0
    END
      ) AS emesis
  FROM
    `physionet-data.eicu_crd.intakeoutput`
    GROUP BY
    patientunitstayid
    )
SELECT
  patientunitstayid
FROM
  blood
UNION DISTINCT
SELECT
  patientunitstayid
FROM
  others
WHERE
  gastric >= 500
  OR stool >= 500
  OR emesis >= 500