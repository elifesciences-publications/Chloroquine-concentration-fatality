;; Chloroquine PK model
;; Based on plasma data from healthy volunteers
;; Pukrittayakamee et al, Antimicrob Agents Chemother, 2014

$PROBLEM    2 ;
$INPUT      ID TIME AMT DV CMT MDV EVID OCC WGT
$DATA       data.csv IGNORE=#
$ABBREVIATED COMRES=2
$SUBROUTINE ADVAN6 TOL=6

$MODEL
COMP=(1)
COMP=(2)
COMP=(3)
COMP=(4)
COMP=(5)
COMP=(6)
COMP=(7)

$PK
TVCL  = THETA(1)*((WGT/70)**0.75)
CL  = TVCL*EXP(ETA(1))

TVV2  = THETA(2)*((WGT/70)**1)
V2  = TVV2* EXP(ETA(2))

TVMT  = THETA(3)
MT  = TVMT* EXP(ETA(3))

TVQ   = THETA(4)*((WGT/70)**0.75)
Q   = TVQ * EXP(ETA(4))

TVV3  = THETA(5)*((WGT/70)**1)
V3  = TVV3* EXP(ETA(5))

TVQ2   = THETA(6)*((WGT/70)**0.75)
Q2   = TVQ2 * EXP(ETA(6))

TVV4  = THETA(7)*((WGT/70)**1)
V4  = TVV4* EXP(ETA(7))

TVF1  = THETA(8)
F1  = TVF1* EXP(ETA(8))

nn=3
KTR = (nn+1)/MT

K15 = KTR
K56 = KTR
K67 = KTR
K72 = KTR
K20 = CL/V2
K23 = Q/V2
K32 = Q/V3
K24 = Q2/V2
K42 = Q2/V4

S2=V2

IF(NEWIND.LE.1) THEN
COM(1)=-1
COM(2)=-1
ENDIF

$DES
DADT(1) = - K15*A(1)
DADT(2) = K72*A(7) + K32*A(3) + K42*A(4) - K23*A(2) - K24*A(2) - K20*A(2)
DADT(3) = K23*A(2) - K32*A(3)
DADT(4) = K24*A(2) - K42*A(4)
DADT(5) = K15*A(1) - K56*A(5)
DADT(6) = K56*A(5) - K67*A(6)
DADT(7) = K67*A(6) - K72*A(7)

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
 41.5  ; 1.CL
 2280  ; 2.V2
 1.26  ; 3.MT
 22.3  ; 4.Q
 7620  ; 5.V3
 191   ; 6.Q2
 3700  ; 7.V4
 1     ; 8.F1

$OMEGA
 0.0236    ; 1.CL
 0.0685    ; 2.V2
 0.151     ; 3.MT
 0 FIX     ; 4.Q
 0 FIX     ; 5.V3
 0 FIX     ; 6.Q2
 0 FIX     ; 7.V4
 0.0111    ; 8.F1

$SIGMA 0.000       ;Simulating without residual error

$SIM (12345) (54321) ONLYSIM SUBPROBLEMS=1000

$TABLE ID TIME AMT CMT IPRED IRES IWRES CWRES MDV EVID OCC WGT CB_CQ1 CB_CQ2 CMAX1 CMAX2
NOPRINT NOTITLE NOLABEL FILE=mytab2_CQ
