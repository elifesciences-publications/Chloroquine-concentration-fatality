;; Chloroquine PK model
;; Based on whole blood data from malaria patients
;; Hoglund et al, Malar J, 2016

$PROBLEM    1 ;
$INPUT      ID TIME AMT DV CMT MDV EVID OCC WGT
$DATA       data.csv IGNORE=#
$ABBREVIATED COMRES=2
$SUBROUTINE ADVAN6 TOL=6

$MODEL
COMP=(1)
COMP=(2)
COMP=(3)
COMP=(4)

$PK
TVCL  = THETA(1)*((WGT/52)**0.75)
CL  = TVCL* EXP(ETA(1))

TVV2  = THETA(2)*((WGT/52)**1)
V2  = TVV2* EXP(ETA(2))

TVMT  = THETA(3)
MT  = TVMT* EXP(ETA(3))

TVQ   = THETA(4)*((WGT/52)**0.75)
Q   = TVQ * EXP(ETA(4))  

TVV3  = THETA(5)*((WGT/52)**1)
V3  = TVV3* EXP(ETA(5))

TVF1  = THETA(6)
F1  = TVF1* EXP(ETA(6))

CLM = 0.176*CL

nn=1
KTR = (nn+1)/MT

K14 = KTR
K42 = KTR
K20 = (CL+CLM)/V2
K23 = Q/V2
K32 = Q/V3

S2=V2

IF(NEWIND.LE.1) THEN
COM(1)=-1
COM(2)=-1
ENDIF

$DES
DADT(1) = - K14*A(1)
DADT(2) = K42*A(4) + K32*A(3) - K23*A(2) - K20*A(2)
DADT(3) = K23*A(2) - K32*A(3)
DADT(4) = K14*A(1) - K42*A(4)

CT=A(2)/V2
IF(TIME.LE.24.AND.CT.GT.COM(1)) THEN
COM(1)=CT
ENDIF

IF(TIME.GT.24.AND.CT.GT.COM(2)) THEN
COM(2)=CT
ENDIF

$ERROR
IPRED  =  A(2)/V2                         ;mg/L
CB_CQ1 = (A(2)/V2)*1000 / 319.872         ;uM
CB_CQ2 = (A(2)/V2)*1000                   ;ng/mL

W     = SQRT(SIGMA(1,1))

IRES  = DV - IPRED
IWRES = IRES/W

Y = IPRED + W*EPS(1)

CMAX1  = COM(1)*1000 / 319.872            ;uM
CMAX2  = COM(2)*1000 / 319.872            ;uM

$THETA
 6.13  ; 1.CL
 468   ; 2.V2
 0.953 ; 3.MT
 37.7  ; 4.Q
 1600  ; 5.V3
 1     ; 6.F1

$OMEGA
 0 FIX  ; 1.CL
 0 FIX  ; 2.V2
 0 FIX  ; 3.MT
 0 FIX  ; 4.Q
 0.0392 ; 5.V3
 0.0369 ; 6.F1

$SIGMA
 0     ;Simulation without residual error

$SIM (12345) (54321) ONLYSIM SUBPROBLEMS=1000

$TABLE ID TIME AMT CMT IPRED IRES IWRES CWRES MDV EVID OCC WGT CB_CQ1 CB_CQ2 CMAX1 CMAX2
NOPRINT NOTITLE NOLABEL FILE=mytab1_CQ
