WITH peakcr 
     AS (SELECT * 
         FROM   (SELECT patientunitstayid, 
                        labresultoffset AS peakcreat7d_offset, 
                        labresult AS peakcreat7d, 
                        Row_number() 
                          OVER ( 
                            partition BY patientunitstayid 
                            ORDER BY labresult DESC) AS position 
                 FROM   `physionet-data.eicu_crd.lab` 
                 WHERE  labname LIKE 'creatinine%' 
                        AND labresultoffset >= 0 
                        AND labresultoffset <= 10080 
                 GROUP  BY patientunitstayid, 
                           labresultoffset, 
                           labresult 
                 ORDER  BY patientunitstayid, 
                           labresultoffset) AS temp 
         WHERE  position = 1) 
SELECT patientunitstayid, 
       peakcreat7d, 
       peakcreat7d_offset, 
       ( unitdischargeoffset - peakcreat7d_offset ) AS 
       peakcreat7d_to_discharge_offsetgap 
FROM   `physionet-data.eicu_crd.patient`
       LEFT OUTER JOIN peakcr USING (patientunitstayid)
ORDER  BY patientunitstayid 


