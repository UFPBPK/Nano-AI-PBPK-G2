NanoPBPK.code <- '

$PARAM @annotated

// Cardiac Output and Blood flow
QCC            : 16.5   :L/h/kg^0.75,                    Cardio output,                  (Brown, 1997)
QLC            : 0.02   :unitless,                       Fraction blood flow to liver,   (Brown, 1997, Table 23)
QLuC           : 1      :unitless,                       Fraction blood flow to lung,    (Brown, 1997, Table 23)
QKC            : 0.091  :unitless,                       Fraction blood flow to kidney,  (Brown, 1997, Table 23)
QBrC           : 0.033  :unitless,                       Fraction blood flow to brain,   (Brown, 1997, Table 23)
QSC            : 0.011  :unitless,                       Fraction blood flow to spleen,  (Lin, 2008; Davies and Morries, 1993)
QMC            : 0.159  :unitless,                       Fraction blood flow to muschle, (Brown, 1997, Table 23)
QTC            : 0.01  :unitless,                       Fraction blood flow to tumor,   

// Tissue Volume
BW             : 0.02   :kg,                             Body weight                    
VLC            : 0.055  :unitless,                       Fraction liver tissue,          (Brown, 1997, Table 21)
VLuC           : 0.007  :unitless,                       Fraction lung tissue,           (Brown, 1997, Table 21)
VKC            : 0.017  :unitless,                       Fraction kidney tissue,         (Brown, 1997, Table 21)
VBrC           : 0.017  :unitless,                       Fraction brain tissue,          (Brown, 1997, Table 21)
VSC            : 0.005  :unitless,                       Fraction spleen tissue,         (Lin, 2008; Davies and Morries, 1993)
VBldC          : 0.06   :unitless,                       Fraction blood,                 (Chen, 2015)
VPlasC         : 0.0355 :unitless,                       Fraction plasma,                (Davies and Morris, 1993
VMC            : 0.384  :unitless,                       Fraction muscle tissue,         (Brown, 1997, Table 21)
VTC            : 0.024   :unitless,                       Fraction tumor,                 (Sykes et al., 2014; Wilhelm et al., 2016)


//Blood volume fraction in organs and tissues
BVL            : 0.31   :unitless,                       Liver,                          (Brown, 1997, Table 30)
BVBr           : 0.03   :unitless,                       Brain,                          (Brown, 1997, Table 30)
BVK            : 0.24   :unitless,                       Kidney,                         (Brown, 1997, Table 30)
BVS            : 0.17   :unitless,                       Spleen,                         (Brown, 1997, Table 30)
BVLu           : 0.5    :unitless,                       lungs,                          (Brown, 1997, Table 30)
BVM            : 0.04   :unitless,                       muscle,                         (Brown, 1997, Table 30)
BVR            : 0.04   :unitless,                       Rest of body (assumed to equal to muscle), (Brown, 1997, Table 30)
BVT            : 0.01   :unitless,                       tumor,                          fitted

//Partition coefficients(PC, tissue:plasma)
PL             : 0.08   :unitless,                       liver,                          (Lin, 2016)
PK             : 0.15   :unitless,                       kidney,                         (Lin, 2016)
PBr            : 0.15   :unitless,                       brain,                          (Lin, 2016)
PS             : 0.15   :unitless,                       spleen,                         (Lin, 2016)
PLu            : 0.15   :unitless,                       lungs,                          (Lin, 2016)
PH             : 0.15   :unitless,                       heart,                          (Lin, 2016)
PM             : 0.15   :unitless,                       muscle,                         (Lin, 2016)
PR             : 0.15   :unitless,                       rest of body,                   (Lin, 2016)
PT             : 0.00377   :unitless,                       tumor,                          fitted

//Membrane-limited permeability coefficient constants
PALC           : 0.001  :unitless,                       liver,                          (Lin, 2016)
PABrC          : 0.000001 :unitless,                     brain,                          (Lin, 2016)
PAKC           : 0.01   :unitless,                       kidney,                         (Lin, 2016)
PASC           : 0.15   :unitless,                       spleen,                         (Lin, 2016)
PALuC          : 0.001  :unitless,                       lung,                           (Lin, 2016)
PAMC           : 0.00005:unitless,                       muscle,                         (Lin, 2016)
PARC           : 0.00005:unitless,                       rest of body,                   (Lin, 2016)
PATC           : 0.00377  :unitless,                       tumor,                          fitted

//Endocytic parameters; RES represent endocytic/phagocytic cells
KLRESrelease   : 0.0015   :1/h,                liver,                          Release rate constant of phyagocytic cells 
KLRESmax       : 0.3      :1/h,                liver,                          Maximmum uptake rate constant of phyagocytic cells
KLRES50        : 48       :h,                  liver,                          Time reaching half maximum uptake rate
KLRESn         : 5        :unitless,           liver,                          Hill coefficient

KSRESrelease   : 0.001    :1/h,                spleen,                         Release rate constant of phyagocytic cells 
KSRESmax       : 5        :1/h,                spleen,                         Maximmum uptake rate constant of phyagocytic cells
KSRES50        : 36       :h,                  spleen,                         Time reaching half maximum uptake rate
KSRESn         : 5        :unitless,           spleen,                         Hill coefficient

KKRESrelease   : 0.001    :1/h,                kidney,                         Release rate constant of phyagocytic cells 
KKRESmax       : 0.12     :1/h,                kidney,                         Maximmum uptake rate constant of phyagocytic cells
KKRES50        : 48       :h,                  kidney,                         Time reaching half maximum uptake rate
KKRESn         : 5        :unitless,           kidney,                         Hill coefficient

KLuRESrelease  : 0.003    :1/h,               lung,                           Release rate constant of phyagocytic cells 
KLuRESmax      : 0.085    :1/h,               lung,                           Maximmum uptake rate constant of phyagocytic cells
KLuRES50       : 48       :h,                 lung,                           Time reaching half maximum uptake rate
KLuRESn        : 5        :unitless,          lung,                           Hill coefficient

KMRESrelease   : 0.005    :1/h,                muscle,                         Release rate constant of phyagocytic cells 
KMRESmax       : 0.4      :1/h,                muscle,                         Maximmum uptake rate constant of phyagocytic cells
KMRES50        : 48       :h,                  muscle,                         Time reaching half maximum uptake rate
KMRESn         : 5        :unitless,           muscle,                         Hill coefficient

KRRESrelease   : 0.005    :1/h,                rest of body,                   Release rate constant of phyagocytic cells 
KRRESmax       : 0.4      :1/h,                rest of body,                   Maximmum uptake rate constant of phyagocytic cells
KRRES50        : 48       :h,                  rest of body,                   Time reaching half maximum uptake rate
KRRESn         : 5        :unitless,           rest of body,                   Hill coefficient

KTRESrelease   : 0.0071    :1/h,                tumor, Release rate constant of phyagocytic cells 
KTRESmax       : 0.00198      :1/h,                tumor, Maximmum uptake rate constant of phyagocytic cells
KTRES50        : 0.00002       :h,                  tumor, Time reaching half maximum uptake rate
KTRESn         : 0.61       :unitless,           tumor, Hill coefficient

//Excretion parameters
KbileC         :  0.00003  :L/hr/kg^0.75,       Bile clearance
KurineC        :  0.000003 :L/hr/kg^0.75,       Urine clearance

$MAIN
double QRC     = 1 - (QLC + QKC  + QSC + QTC + QBrC + QMC);                          //Fraction of blood flow to rest of body
double VRC     = 1 - (VLC + VLuC + VKC  + VSC + VBrC + VMC + VTC + VPlasC);                          //Tissue volume of rest of body

double QC      = QCC * pow(BW, 0.75);                                                    //L/h, Cardiac output (adjusted for plasma)
double QL      = QC * QLC;                                                           //L/h, Blood flow to liver
double QBr     = QC * QBrC;                                                         //L/h, Blood flow to brain
double QK      = QC * QKC;                                                           //L/h, Blood flow to kidney
double QS      = QC * QSC;                                                          //L/h, Blood flow to spleen
double QM      = QC * QMC;                                                           //L/h, Blood flow to muscle
double QT      = QC * QTC;                                                           //L/h, Blood flow to tumor
double QR      = QC * QRC;                                                           //L/h, Blood flow to the rest of body

//Tissue volumes
double VL      = BW * VLC;                                                           //L, Liver volume
double VBr     = BW * VBrC;                                                          //L, Brain volume
double VK      = BW * VKC;                                                           //L, Kidney volume
double VM      = BW * VMC;
double VS      = BW * VSC;                                                          //L, spleen volume
double VLu     = BW * VLuC;                                                          //L, lung volume
double VR      = BW * VRC;                                                           //L, Rest of body
double VT      = BW * VTC;                                                           //L, Tumor
double VBlood  = BW * VBldC;                                                      //L, Blood
double VPlasma = BW * VPlasC;                                                    //L, Plasma

double VLb     = VL * BVL; 							                                            //Weight/volume of capillary blood in liver compartment
double VLt     = VL - VLb; 							                                            //Weight/volume of tissue in liver compartment
double VBrb    = VBr * BVBr; 						                                            //Weight/volume of capillary blood in brain compartment
double VBrt    = VBr - VBrb; 						                                            //Weight/volume of tissue in brain compartment
double VKb     = VK * BVK; 							                                            //Weight/volume of capillary blood in kidney compartment
double VKt     = VK - VKb; 							                                            //Weight/volume of tissue in kidney compartment
double VSb     = VS * BVS; 							                                          //Weight/volume of capillary blood in spleen compartment
double VSt     = VS - VSb; 							                                          //Weight/volume of tissue in spleen compartment
double VLub    = VLu * BVLu; 						                                            //Weight/volume of capillary blood in lung compartment
double VLut    = VLu - VLub; 						                                            //Weight/volume of tissue in lung compartment
double VMb     = VM * BVM; 						                                              //Weight/volume of capillary blood in muscle compartment
double VMt     = VM - VMb; 							                                            //Weight/volume of tissue in muscle compartment
double VRb     = VR * BVR; 							                                            //Weight/volume of capillary blood in rest of body compartment
double VRt     = VR - VRb; 							                                            //Weight/volume of tissue in rest of body compartment
double VTb     = VT * BVT; 							                                            //Weight/volume of capillary blood in tumor compartment
double VTt     = VT - VTb; 							                                            //Weight/volume of tissue in tumor compartment

//Permeability coefficient-surface area cross-product (L/h)
double PAL     = PALC  * QL; 						                                            //L/h, Liver
double PABr    = PABrC * QBr; 						                                          //L/h, Brain
double PAK     = PAKC  * QK; 						                                            //L/h, Kidneys
double PAS     = PASC  * QS; 						                                          //L/h, Spleen
double PALu    = PALuC * QC;  						                                          //L/h, Lungs
double PAM     = PAMC  * QM; 						                                            //L/h, Muscle
double PAR     = PARC  * QR; 						                                            //L/h, Rest of body
double PAT     = PATC  * QT; 						                                            //L/h, Tumor

//Endocytosis rate (1/h)
double KLRESUP = ((KLRESmax)*pow(TIME, KLRESn))/(pow(KLRES50, KLRESn) + pow(TIME, KLRESn)); 	    
double KSRESUP = ((KSRESmax)*pow(TIME, KSRESn))/(pow(KSRES50, KSRESn) + pow(TIME, KSRESn)); 	    //Lungs
double KKRESUP = ((KKRESmax)*pow(TIME, KKRESn))/(pow(KKRES50, KKRESn) + pow(TIME, KKRESn)); 	    //Lungs
double KLuRESUP = ((KLuRESmax)*pow(TIME, KLuRESn))/(pow(KLuRES50, KLuRESn) + pow(TIME, KLuRESn)); 	    //Lungs
double KMRESUP  = ((KMRESmax)*pow(TIME, KMRESn))/(pow(KMRES50, KMRESn) + pow(TIME, KMRESn));	          //Muscle
double KRRESUP  = ((KRRESmax)*pow(TIME, KRRESn))/(pow(KRRES50, KRRESn) + pow(TIME, KRRESn));		        //Rest of body

double KTRESUP1 = (KTRESmax*pow(TIME,KTRESn))/(pow(KTRES50, KTRESn) + pow(TIME, KTRESn));		              
double Kbile   = KbileC * pow(BW, 0.75);               
double Kurine  = KurineC * pow(BW, 0.75);

$CMT AA AV ALub ALut ALuRES ABrb ABrt AMb AMt AMRES ARb ARt ARRES 
AKb AKt AKRES Aurine ASb ASt ASRES ALb ALt ALRES Abile ATb ATt ATRES ATRESUP ATRESrel ADOSE 
AUCTumor 

$ODE

double APlasma  = AA + AV;
double ABlood   = AA + AV;
double ALung    = ALub+ALut+ALuRES;
double ALungt   = ALut+ALuRES;
double Arestall = ARb+ARt+ARRES;
double Aresttissue = ARt+ARRES;
double AKidney  = AKb+AKt+AKRES;
double AKidneyt = AKt+AKRES;
double ABrain   = ABrb + ABrt;
double ASpleen  = ASb+ASt+ASRES;
double ASpleent = ASt+ASRES;
double AMuscle  = AMb+AMt+AMRES;
double AMusclet = AMt+AMRES;
double ALiver   = ALb+ALt+ALRES;
double ALivert  = ALt+ALRES;
double ATumor   = ATb+ATt+ATRES;
double ATumort  = ATt+ATRES;


double CA       = AA/(VPlasma * 0.2);
double CV       = AV / (VPlasma * 0.8);
double CPlasma  = APlasma/VPlasma;
double CBlood   = ABlood/VBlood;
double CVLu     = ALub / VLub;
double CVL      = ALb/VLb;
double CLut     = ALut / VLut;
double CLung    = (ALub + ALut + ALuRES)/VLu;
double CLungt   = (ALut+ALuRES)/VLut;
double CVBr     = ABrb/VBrb;
double CBrt     = ABrt/VBrt;
double CBrain   = ABrain / VBr;
double CVK      = AKb/VKb;
double CVM      = AMb/VMb;
double CMt      = AMt/VMt;
double CMuscle  = (AMb+AMt+AMRES)/VM;
double CMusclet = (AMt+AMRES)/VMt;
double CVR      = ARb/VRb;
double CRt      = ARt/VRt;
double Crestall = (ARb+ARt+ARRES)/VR;
double Cresttissue = (ARt+ARRES)/VRt;
double CKt      = AKt/VKt;
double CKidney  = (AKb+AKt+AKRES)/VK;
double CKidneyt = (AKt+AKRES)/VKt;
double CVS      = ASb/VSb;
double CSt      = ASt/VSt;
double CSpleen  = (ASb+ASt+ASRES)/VS;
double CSpleent = (ASt+ASRES)/VSt;
double CLt      = ALt/VLt;
double CLiver = (ALb+ALt+ALRES)/VL;
double CLivert = (ALt+ALRES)/VLt;
double CVT     = ATb/VTb;
double CTt      = ATt/VTt;
double CTumort  = (ATt+ATRES)/VTt;
double CTumor   = (ATb+ATt+ATRES)/VT;


double RA      = QC * CVLu - QC * CA; 
double RV      = QL * CVL + QBr  *CVBr + QK * CVK + QM * CVM + QR * CVR + QT * CVT - QC * CV;
double RLub    =  QC * (CV - CVLu) - PALu * CVLu + (PALu * CLut)/ PLu; 
double RLut    = PALu * CVLu - (PALu * CLut)/ PLu - KLuRESUP * ALut + KLuRESrelease * ALuRES;
double RLuRES  = KLuRESUP * ALut - KLuRESrelease * ALuRES;
double RBrb    = QBr *(CA - CVBr) - PABr * CVBr + (PABr * CBrt)/ PBr;
double RBrt    = PABr * CVBr - (PABr * CBrt)/ PBr;
double RMb      = QM*(CA-CVM) - PAM*CVM + (PAM*CMt)/PM;
double RMt      = PAM*CVM - (PAM*CMt)/PM - KMRESUP*AMt + KMRESrelease*AMRES;
double RMRES    = KMRESUP*AMt-KMRESrelease*AMRES;
double RRb      = QR*(CA-CVR) - PAR*CVR + (PAR*CRt)/PR;
double RRt      = PAR*CVR - (PAR*CRt)/PR - KRRESUP*ARt + KRRESrelease*ARRES;
double RRRES    = KRRESUP*ARt-KRRESrelease*ARRES;
double Rurine   = Kurine*CVK;
double RKb      = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK - Rurine;
double RKt      = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES; 
double RKRES    = KKRESUP*AKt-KKRESrelease*AKRES;
double RSb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS; 
double RSt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES;
double RSRES    = KSRESUP*ASt-KSRESrelease*ASRES;
double Rbile    = Kbile*CLt ; 
double RLb      = QL*(CA-CVL) + QS*CVS - PAL*CVL + (PAL*CLt)/PL - KLRESUP*ALb + KLRESrelease*ALRES;
double RLt      = PAL*CVL - (PAL*CLt)/PL - Rbile;
double RLRES    = KLRESUP*ALb-KLRESrelease*ALRES;
double RTb      = QT*(CA-CVT) - PAT*CVT + (PAT*CTt)/PT; 
double RTt      = PAT*CVT - (PAT*CTt)/PT - KTRESUP1*ATt + KTRESrelease*ATRES;
double RTRES    = KTRESUP1*ATt-KTRESrelease*ATRES;
double RTRESUP  = KTRESUP1*ATt;
double RTRESrel = KTRESrelease*ATRES;

//Concentration of the chemical in blood compartment
//CA: Arterial blood concentration (mg/L or ug/ml)
dxdt_AA     = RA;
dxdt_AV     = RV;
dxdt_ALub   = RLub;
dxdt_ALut   = RLut;
dxdt_ALuRES = RLuRES;
dxdt_ABrb   = RBrb;
dxdt_ABrt   = RBrt;
dxdt_AMb    = RMb;
dxdt_AMt    = RMt;
dxdt_AMRES  = RMRES;
dxdt_ARb    = RRb;
dxdt_ARt    = RRt;
dxdt_ARRES  = RRRES;
dxdt_AKb    = RKb;
dxdt_AKt    = RKt;
dxdt_AKRES  = RKRES;
dxdt_Aurine = Rurine;
dxdt_Abile  = Rbile;
dxdt_ASb    = RSb;
dxdt_ASt    = RSt;
dxdt_ASRES  = RSRES;
dxdt_ALb    = RLb;
dxdt_ALt    = RLt;
dxdt_ALRES  = RLRES;
dxdt_Abile  = Rbile;
dxdt_ATb    = RTb;
dxdt_ATt         = RTt;
dxdt_ATRES       = RTRES;
dxdt_ATRESUP     = RTRESUP;
dxdt_ATRESrel    = RTRESrel;
dxdt_ADOSE       = 0;
dxdt_AUCTumor    = CTumor;

// {Mass balance equations}
double Tmass = APlasma + ALiver + ABrain + AKidney + ALung + Arestall + AMuscle + ASpleen + Abile + Aurine + ATumor;
double Bal   = ADOSE-Tmass;

$TABLE
capture AUCT     = AUCTumor;
capture Tumor    = ATt+ATRES;
capture Lung     = CLung;
capture Liver    = CLiver;
capture Kidney   = CKidney;
capture Spleen   = CSpleen;
capture BAL      = Bal;
'