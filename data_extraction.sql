-- MIMIC-IV (v2.0)
WITH sofa AS
(
  SELECT stay_id
    , starttime, endtime
    , respiration_24hours as respiration
    , coagulation_24hours as coagulation
    , liver_24hours as liver
    , cardiovascular_24hours as cardiovascular
    , cns_24hours as cns
    , renal_24hours as renal
    , sofa_24hours as sofa_score
  FROM `physionet-data.mimiciv_derived.sofa`
  WHERE sofa_24hours >= 2
)
, s1 as
(
  SELECT 
    soi.subject_id
    , soi.stay_id
    -- suspicion columns
    , soi.ab_id
    , soi.antibiotic
    , soi.antibiotic_time
    , soi.culture_time
    , soi.suspected_infection
    , soi.suspected_infection_time
    , soi.specimen
    , soi.positive_culture
    -- sofa columns
    , starttime, endtime
    , respiration, coagulation, liver, cardiovascular, cns, renal
    , sofa_score
    -- All rows have an associated suspicion of infection event
    -- Therefore, Sepsis-3 is defined as SOFA >= 2.
    -- Implicitly, the baseline SOFA score is assumed to be zero, as we do not know
    -- if the patient has preexisting (acute or chronic) organ dysfunction 
    -- before the onset of infection.
    , sofa_score >= 2 and suspected_infection = 1 as sepsis3
    -- subselect to the earliest suspicion/antibiotic/SOFA row
    , ROW_NUMBER() OVER
    (   PARTITION BY soi.stay_id
        ORDER BY suspected_infection_time, antibiotic_time, culture_time, endtime
    ) AS rn_sus
  FROM `physionet-data.mimiciv_derived.suspicion_of_infection` as soi
  INNER JOIN sofa
    ON soi.stay_id = sofa.stay_id 
    AND sofa.endtime >= DATETIME_SUB(soi.suspected_infection_time, INTERVAL 48 HOUR)
    AND sofa.endtime <= DATETIME_ADD(soi.suspected_infection_time, INTERVAL 24 HOUR)
  -- only include in-ICU rows
  WHERE soi.stay_id is not null
),
s2 as
(
    select subject_id, 
    max(case when positive_culture = 1 and specimen in ('BLOOD', 'Blood (EBV)','Blood (LYME)','Blood (Toxo)','BLOOD CULTURE','SEROLOGY/BLOOD',
    'Blood (Malaria)','CATHETER TIP-IV','Immunology (CMV)','BLOOD CULTURE - NEONATE','Stem Cell - Blood Culture','BLOOD CULTURE (POST-MORTEM)',
    'PERIPHERAL BLOOD LYMPHOCYTES','BLOOD CULTURE ( MYCO/F LYTIC BOTTLE)','Blood (CMV AB)') then 1 else 0 end) as blood,
    max(case when positive_culture = 1 and specimen in ('SPUTUM','Mini-BAL','THROAT CULTURE','BRONCHIAL BRUSH','PLEURAL FLUID','THROAT FOR STREP','TRACHEAL ASPIRATE',
    'BRONCHIAL WASHINGS','BRONCHOALVEOLAR LAVAGE','BRONCHIAL BRUSH - PROTECTED','Influenza A/B by DFA - Bronch Lavage','RAPID RESPIRATORY VIRAL ANTIGEN TEST',
    'Rapid Respiratory Viral Screen & Culture') then 1 else 0 end) as resp,
    max(case when positive_culture = 1 and specimen in ('CSF;SPINAL FLUID') then 1 else 0 end) as csf,
    max(case when positive_culture = 1 and specimen in ('BILE', 'STOOL','FECAL SWAB','PERITONEAL FLUID','STOOL (RECEIVED IN TRANSPORT SYSTEM)') then 1 else 0 end) as gi,
    max(case when positive_culture = 1 and specimen in ('URINE','URINE,KIDNEY','ANORECTAL/VAGINAL CULTURE','URINE,SUPRAPUBIC ASPIRATE') then 1 else 0 end) as gu,
    max(case when positive_culture = 1 and specimen in ('EAR', 'EYE', 'SWAB', 'Swab', 'WORM', 'BIOPSY','TISSUE','ABSCESS','Isolate','ASPIRATE','ARTHROPOD',
    'CRE Screen','IMMUNOLOGY','BONE MARROW','FLUID WOUND','FLUID,OTHER','JOINT FLUID','MRSA SCREEN','FOOT CULTURE','FOREIGN BODY',
    'SWAB, R/O GC','VIRAL CULTURE','DIALYSIS FLUID','NAIL SCRAPINGS','SKIN SCRAPINGS','BLOOD BAG FLUID','NEOPLASTIC BLOOD','Staph aureus swab',
    'POSTMORTEM CULTURE','C, E, & A Screening','Touch Prep/Sections','Influenza A/B by DFA','CORNEAL EYE SCRAPINGS','Swab R/O Yeast Screen',
    'PROSTHETIC JOINT FLUID','Infection Control Yeast','SCOTCH TAPE PREP/PADDLE','BONE MARROW - CYTOGENETICS','Foreign Body - Sonication Culture',
    'VIRAL CULTURE: R/O CYTOMEGALOVIRUS','VIRAL CULTURE:R/O HERPES SIMPLEX VIRUS','FLUID RECEIVED IN BLOOD CULTURE BOTTLES',
    'DIRECT ANTIGEN TEST FOR VARICELLA-ZOSTER VIRUS','Direct Antigen Test for Herpes Simplex Virus Types 1 & 2') then 1 else 0 end) as other  
  from s1
  group by subject_id
),
cohort as
(SELECT 
subject_id, stay_id
-- note: there may be more than one antibiotic given at this time
, antibiotic_time, antibiotic
-- culture times may be dates, rather than times
, culture_time, specimen
, suspected_infection_time
-- endtime is latest time at which the SOFA score is valid
, endtime as sofa_time
, sofa_score
, respiration, coagulation, liver, cardiovascular, cns, renal
, sepsis3
, case when suspected_infection_time <= endtime then endtime else suspected_infection_time end as sepsis_time
FROM s1
WHERE rn_sus = 1
and sepsis3 = true
),
uo as
(
  select i.stay_id, sum(u.urineoutput) as uo_firstday
  from `physionet-data.mimiciv_derived.icustay_detail` i
  left join  `physionet-data.mimiciv_derived.urine_output` u
  on u.stay_id = i.stay_id
  where u.charttime >= i.icu_intime
  and u.charttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.stay_id
  order by i.stay_id
),
gcs as
(
  select i.stay_id, min(g.gcs) as gcs_min
  from  `physionet-data.mimiciv_derived.icustay_detail` i
  left join `physionet-data.mimiciv_derived.gcs` g
  on g.stay_id = i.stay_id
  where g.charttime >= i.icu_intime
  and g.charttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.stay_id
  order by i.stay_id
),
rrt as
(
  select i.stay_id, 
  max(case when r.dialysis_type in ('Peritoneal','IHD','SCUF','CVVH','CVVHD','CVVHDF') then 1 else 0 end) as rrt_flag,
  MAX(case when r.dialysis_type = 'Peritoneal' then 1 else 0 end) as peritoneal,
  MAX(case when r.dialysis_type = 'IHD' then 1 else 0 end) as ihd,
  MAX(case when r.dialysis_type = 'SCUF' then 1 else 0 end) as scuf,
  MAX(case when r.dialysis_type = 'CVVH' then 1 else 0 end) as cvvh,
  MAX(case when r.dialysis_type = 'CVVHD' then 1 else 0 end) as cvvhd,
  MAX(case when r.dialysis_type = 'CVVHDF' then 1 else 0 end) as cvvhdf
  from `physionet-data.mimiciv_derived.icustay_detail` i 
  left join `physionet-data.mimiciv_derived.rrt` r
  on r.stay_id = i.stay_id
  where r.charttime >= i.icu_intime
  and r.charttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.stay_id
  order by i.stay_id
),
vent as
(
  select i.stay_id, 
  MAX(case when v.ventilation_status in ('Tracheostomy','InvasiveVent') then 1 else 0 end) as invasive,
  MAX(case when v.ventilation_status = 'NonInvasiveVent' then 1 else 0 end) as noninvasive,
  MAX(case when v.ventilation_status = 'SupplementalOxygen' then 1 else 0 end) as ox,
  MAX(case when v.ventilation_status = 'HFNC' then 1 else 0 end) as hfhc
  from `physionet-data.mimiciv_derived.icustay_detail` i
  left join  `physionet-data.mimiciv_derived.ventilation` v
  on v.stay_id = i.stay_id
  where v.starttime >= i.icu_intime
  and v.starttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.stay_id
  order by i.stay_id
),
dobutamine AS
(SELECT i.stay_id, 
MAX(CASE WHEN db.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
and db.vaso_rate is not null THEN 1 ELSE 0 END) AS dobutamine_flag,
MAX(CASE WHEN db.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN db.vaso_rate ELSE NULL END) AS max_dobutamine
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.dobutamine` db
ON db.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
dopamine AS
(SELECT i.stay_id, 
MAX(CASE WHEN d.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
and d.vaso_rate is not null THEN 1 ELSE 0 END) AS dopamine_flag,
MAX(CASE WHEN d.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN d.vaso_rate ELSE NULL END) AS max_dopamine
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.dopamine` d
ON d.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
epinephrine AS
(SELECT i.stay_id, 
MAX(CASE WHEN e.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
and e.vaso_rate is not null THEN 1 ELSE 0 END) AS epinephrine_flag,
MAX(CASE WHEN e.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN e.vaso_rate ELSE NULL END) AS max_epinephrine
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.epinephrine` e
ON e.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
norepinephrine AS
(SELECT i.stay_id, 
MAX(CASE WHEN e.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
and e.vaso_rate is not null THEN 1 ELSE 0 END) AS norepinephrine_flag,
MAX(CASE WHEN e.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN e.vaso_rate ELSE NULL END) AS max_norepinephrine
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.norepinephrine` e
ON e.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
phenylephrine AS
(SELECT i.stay_id, 
MAX(CASE WHEN p.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
and p.vaso_rate is not null THEN 1 ELSE 0 END) AS phenylephrine_flag,
MAX(CASE WHEN p.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN p.vaso_rate ELSE NULL END) AS max_phenylephrine
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.phenylephrine` p
ON p.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
vasopressin AS
(SELECT i.stay_id, 
MAX(CASE WHEN v.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
and v.vaso_rate is not null THEN 1 ELSE 0 END) AS vasopressin_flag,
MAX(CASE WHEN v.starttime BETWEEN i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR) 
THEN v.vaso_rate ELSE NULL END) AS max_vasopressin
FROM `physionet-data.mimiciv_derived.icustay_detail` i
LEFT JOIN  `physionet-data.mimiciv_derived.vasopressin` v
ON v.stay_id = i.stay_id
GROUP BY i.stay_id
order by i.stay_id
),
crp as
(
  select i.hadm_id, max(f.crp) as crp_max, min(f.crp) as crp_min
  from `physionet-data.mimiciv_derived.icustay_detail` i
  left join `physionet-data.mimiciv_derived.inflammation` f
  on i.hadm_id = f.hadm_id
  where f.charttime between i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.hadm_id
  order by i.hadm_id   
),
br as
(
  select i.hadm_id, max(b.mch) as mch_max, min(b.mch) as mch_min, max(b.mchc) as mchc_max, min(b.mchc) as mchc_min,
  max(b.mcv) as mcv_max, min(b.mcv) as mcv_min, max(b.rbc) as rbc_max, min(b.rbc) as rbc_min, 
  max(b.rdw) as rdw_max, min(b.rdw) as rdw_min
  FROM  `physionet-data.mimiciv_derived.icustay_detail` i
  left join `physionet-data.mimiciv_derived.complete_blood_count`  b
  on i.hadm_id = b.hadm_id
  where b.charttime between i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.hadm_id
  order by i.hadm_id
),
tn as
(
  select i.hadm_id, 
  max(cm.ck_mb) as ck_max, min(cm.ck_mb) as ck_min, 
  max(cm.ntprobnp) as probnp_max, min(cm.ntprobnp) as probnp_min
  from `physionet-data.mimiciv_derived.icustay_detail` i
  left join `physionet-data.mimiciv_derived.cardiac_marker` cm
  on i.hadm_id = cm.hadm_id
  where cm.charttime between i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
  group by i.hadm_id
  order by i.hadm_id
),
nlr as
(
SELECT
  n.hadm_id, max(n.neutrophils) as neut_max, min(n.neutrophils) as neut_min, max(n.lymphocytes) as lymp_max,
  min(n.lymphocytes) as lymp_min, max(n.nlr) as nlr_max, min(n.nlr) as nlr_min
FROM
(  SELECT hadm_id, charttime, neutrophils, lymphocytes, (neutrophils/lymphocytes) AS nlr
   FROM `physionet-data.mimiciv_derived.blood_differential`
   where neutrophils > 0 AND lymphocytes > 0 
) n
left join `physionet-data.mimiciv_derived.icustay_detail` i
on n.hadm_id = i.hadm_id
where n.charttime between i.icu_intime and DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
group by n.hadm_id
order by n.hadm_id
),
comorbidities as (
  select hadm_id, 
  max(case when diabetes_without_cc = 1 or diabetes_with_cc = 1 then 1 else 0 end) as diabetes,
  max(case when malignant_cancer = 1 or metastatic_solid_tumor = 1 then 1 else 0 end) as cancer
  from `physionet-data.mimiciv_derived.charlson` 
  group by hadm_id
),
nu0 as (
    select i.stay_id, MAX(case when g.itemid = 227079 then 1 else 0 end) as PN_flag,
    MAX(case when g.itemid in (226221,227080) then 1 else 0 end) as EN_flag
    from `physionet-data.mimiciv_derived.icustay_detail` i
    left join `physionet-data.mimiciv_icu.ingredientevents` g
    on i.stay_id = g.stay_id
    where g.starttime >= i.icu_intime and g.starttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
    group by i.stay_id
    order by i.stay_id
),
nud1 as (
    select i.stay_id, 
    SUM(case when g.itemid in (226060, 220412) then g.amount else 0 end) as calories_sum,
    SUM(case when g.itemid = 220454 then g.amount else 0 end) as protein_sum,
    SUM(case when g.itemid = 220490 then g.amount else 0 end) as volume_sum
    from `physionet-data.mimiciv_derived.icustay_detail` i 
    left join `physionet-data.mimiciv_icu.ingredientevents` g 
    on i.stay_id = g.stay_id
    where g.starttime >= i.icu_intime and g.starttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
    group by i.stay_id
    order by i.stay_id
),
citrate1 as (
    select i.stay_id, 
    SUM(case when p1.itemid = 227526 then p1.amount else 0 end)*0.59 as d1_calories_citrate
    from `physionet-data.mimiciv_derived.icustay_detail` i
    left join `physionet-data.mimiciv_icu.inputevents` p1
    on i.stay_id = p1.stay_id
    and p1.starttime >= i.icu_intime and p1.starttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)   
    group by i.stay_id
    order by i.stay_id
),
insulin1 as
(
    select i.stay_id,
    sum(case when p.itemid in (223257,223258,223259,223260,223261,223262,229299,229619) then p.amount else 0 end) as d1_insulin
    from `physionet-data.mimiciv_derived.icustay_detail` i
    left join `physionet-data.mimiciv_icu.inputevents` p
    on i.stay_id = p.stay_id
    where p.starttime >= i.icu_intime and p.starttime <= DATETIME_ADD(i.icu_intime, INTERVAL 24 HOUR)
    group by i.stay_id
    order by i.stay_id
),
final as
(select 
i.subject_id, i. hadm_id, i. stay_id, i.gender, i.admission_age, i.race, 
i.dod, i.hospital_expire_flag, i.los_hospital, i.hospstay_seq, 
i.icustay_seq, 
round(datetime_diff(i.icu_outtime, i.icu_intime, minute)/60/24,2) as los_icu, 
round(datetime_diff(i.dischtime, i.icu_intime, minute)/60/24,2) as los_hosp_icu,
round(datetime_diff(i.dod, i.icu_intime, minute)/60/24,2) as los_mortality,
round(datetime_diff(c.antibiotic_time, i.icu_intime, minute)/60,2) as anti_gap,
round(datetime_diff(c.culture_time, i.icu_intime, minute)/60,2) as cul_gap,
round(datetime_diff(c.sepsis_time, i.icu_intime, minute)/60,2) as sepsis_gap,
h.height, w.weight_admit AS weight, icu.first_careunit,
a.admission_type, a.admission_location, a.insurance, a.marital_status,  p.anchor_year_group,
-- comorbidities
co.myocardial_infarct as MI, co.congestive_heart_failure as CHF, co.cerebrovascular_disease as CVD,
co.chronic_pulmonary_disease as cpd, co.renal_disease as ckd,
co1.diabetes, co1.cancer,
-- vital signs
v.heart_rate_min, v.heart_rate_max, v.heart_rate_mean, v.sbp_min, v.sbp_max, v.sbp_mean, v.dbp_min, v.dbp_max, v.dbp_mean,
v.mbp_min, v.mbp_max, v.mbp_mean, v.resp_rate_min, v.resp_rate_max, v.resp_rate_mean, v.temperature_min, v.temperature_max,
v.temperature_mean, v.spo2_min, v.spo2_max, v.spo2_mean,
-- laboratory
-- blood routine
l.wbc_min, l.wbc_max, l.hematocrit_min, l.hematocrit_max, l.hemoglobin_min, l.hemoglobin_max, l.platelets_min, l.platelets_max,
crp.crp_max, crp.crp_min,
br.mch_max, br.mch_min, br.mchc_max, br.mchc_min, br.mcv_max, br.mcv_min, br.rbc_max, 
br.rbc_min, br.rdw_max, br.rdw_min, 
nlr.neut_max, nlr.neut_min, nlr.lymp_max, nlr.lymp_min, nlr.nlr_max, nlr.nlr_min,
-- kidney
l.bun_min, l.bun_max, l.creatinine_min, l.creatinine_max, 
l.bicarbonate_min, l.bicarbonate_max, l.sodium_min, l.sodium_max,
l.potassium_min, l.potassium_max, l.chloride_min, l.chloride_max, 
l.calcium_min, l.calcium_max, l.aniongap_min, l.aniongap_max,
-- liver
l.albumin_min, l.albumin_max, l.globulin_min, l.globulin_max, 
l.total_protein_min, l.total_protein_max, 
l.alt_min, l.alt_max, l.alp_min, l.alp_max, l.ast_min, l.ast_max, 
l.bilirubin_total_min, l.bilirubin_total_max, l.bilirubin_direct_min, l.bilirubin_direct_max, l.bilirubin_indirect_min, l.bilirubin_indirect_max, 
l.ggt_min, l.ggt_max, l.amylase_min, l.amylase_max,
l.ld_ldh_min, l.ld_ldh_max, l.glucose_min, l.glucose_max, 
-- heart
tn.ck_max, tn.ck_min, tn.probnp_max, tn.probnp_min,
-- cogualation
l.pt_min, l.pt_max, l.ptt_min, l.ptt_max, l.inr_min, l.inr_max, 
l.fibrinogen_min, l.fibrinogen_max, l.thrombin_min, l.thrombin_max,
l.d_dimer_min, l.d_dimer_max,
-- blood gas
bg.lactate_min, bg.lactate_max, 
bg.ph_min, bg.ph_max, bg.baseexcess_min as BE_min, bg.baseexcess_max as BE_max, bg.totalco2_min, bg.totalco2_max,
bg1.po2_min, bg1.po2_max, bg1.pco2_min, bg1.pco2_max,
bg1.aado2_min, bg1.aado2_max, bg1.pao2fio2ratio_min as OI_min, bg1.pao2fio2ratio_max as OI_max,
-- uo
uo.uo_firstday,
-- gcs
gcs.gcs_min,
-- treatment
-- rrt
rrt.rrt_flag, rrt.peritoneal, rrt.ihd, rrt.scuf, rrt.cvvh, rrt.cvvhd, rrt.cvvhdf,
-- vent
vent.invasive,vent.noninvasive,vent.ox,vent.hfhc,
-- vasopressor 
dobutamine.dobutamine_flag, dobutamine.max_dobutamine, 
dopamine.dopamine_flag, dopamine.max_dopamine,
epinephrine.epinephrine_flag, epinephrine.max_epinephrine,
norepinephrine.norepinephrine_flag, norepinephrine.max_norepinephrine,
phenylephrine.phenylephrine_flag, phenylephrine.max_phenylephrine,
vasopressin.vasopressin_flag, vasopressin.max_vasopressin,
-- nutrition
nu0.PN_flag, nu0.EN_flag,
(nud1.calories_sum + citrate1.d1_calories_citrate) as d1_calories,
nud1.protein_sum, nud1.volume_sum,
insulin1.d1_insulin,
-- infection
s2.blood, s2.resp, s2.csf, s2.gi, s2.gu, s2.other, 
-- severity score
c.sofa_score, lods.LODS, sapsii.sapsii as SAPSII, sapsiii.apsiii as SAPSIII, oasis.oasis as OASIS,  
ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.hospstay_seq , i.icustay_seq) as seqnum
from `physionet-data.mimiciv_derived.icustay_detail` i
inner join cohort c
on i.stay_id = c.stay_id
left join `physionet-data.mimiciv_derived.height` h 
on i.stay_id = h.stay_id
left join `physionet-data.mimiciv_derived.first_day_weight` w 
on i.stay_id = w.stay_id
left join `physionet-data.mimiciv_icu.icustays` icu
on i.stay_id = icu.stay_id
left join `physionet-data.mimiciv_hosp.admissions` a 
on i.hadm_id = a.hadm_id
left join `physionet-data.mimiciv_hosp.patients` p 
on i.subject_id = p.subject_id
left join `physionet-data.mimiciv_derived.charlson` co 
on i.hadm_id = co.hadm_id 
left join comorbidities co1 
on i.hadm_id = co1.hadm_id
left join `physionet-data.mimiciv_derived.first_day_vitalsign` v
on i.stay_id = v.stay_id
left join `physionet-data.mimiciv_derived.first_day_lab` l
on i.stay_id = l.stay_id
left join `physionet-data.mimiciv_derived.first_day_bg` bg
on i.stay_id = bg.stay_id
left join `physionet-data.mimiciv_derived.first_day_bg_art` bg1
on i.stay_id = bg1.stay_id
left join crp
on i.hadm_id = crp.hadm_id
left join br
on i.hadm_id = br.hadm_id
left join tn
on i.hadm_id = tn.hadm_id
left join nlr
on i.hadm_id = nlr.hadm_id
left join uo
on i.stay_id = uo.stay_id
left join gcs
on i.stay_id = gcs.stay_id
left join rrt
on i.stay_id = rrt.stay_id
left join vent
on i.stay_id = vent.stay_id
left join dobutamine
on i.stay_id = dobutamine.stay_id
left join dopamine
on i.stay_id = dopamine.stay_id
left join epinephrine
on i.stay_id = epinephrine.stay_id
left join norepinephrine
on i.stay_id = norepinephrine.stay_id
left join phenylephrine
on i.stay_id = phenylephrine.stay_id
left join vasopressin
on i.stay_id = vasopressin.stay_id
left join nu0
on i.stay_id = nu0.stay_id
left join nud1
on i.stay_id = nud1.stay_id
left join citrate1
on i.stay_id = citrate1.stay_id
left join insulin1
on i.stay_id = insulin1.stay_id
left join s2
on i.subject_id = s2.subject_id
left join `physionet-data.mimiciv_derived.lods` lods
on i.stay_id = lods.stay_id
left join `physionet-data.mimiciv_derived.sapsii`  sapsii
on i.stay_id = sapsii.stay_id
left join `physionet-data.mimiciv_derived.apsiii` sapsiii
on i.stay_id = sapsiii.stay_id
left join `physionet-data.mimiciv_derived.oasis` oasis
on i.stay_id = oasis.stay_id
)
select final.*
from final
where final.seqnum = 1
and final.admission_age >= 18
and final.los_icu >= 1
order by final.subject_id, final.icustay_seq;
