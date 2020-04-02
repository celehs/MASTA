rm(list = ls())
####### source PP_FPCA function
####### Function for loading/installing necessary package
pkgTest <- function(x){
  if (!require(x,character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest('foreign')

mac = TRUE;type = "lung";
if(mac){
  wkdir  = "~/Google Drive/CanCORS & CRN/OrigData/"
  wkdir2 = "~/Google Drive/CanCORS & CRN/CleanData/"
} else{
  wkdir  = "C:/Users/llian/Google Drive/CanCORS & CRN/OrigData/" 
  wkdir2 = "C:/Users/llian/Google Drive/CanCORS & CRN/CleanData/" 
}

cancors   = read.dta(paste0(wkdir,"data_cancors_ph2.dta"))
crn       = read.dta(paste0(wkdir,"data_crn_ph2.dta"))
crn$case  = sapply(crn$case,function(x) substr(x,4,nchar(x)))
## unlabeled
SEER_med  = read.dta(paste0(wkdir,"d05_recmortality_algorithminput_ver02.dta")) # secondary_mos missing
names(SEER_med)[6] = "analysisfu"

### extract only crc/lung cancer patients
cancors   = cancors[cancors$cancertype==type,]
crn       = crn[crn$cancertype==type,]
SEER_med  = SEER_med[SEER_med$cancertype==type,] # no overlap between two data sets

### codes (nanyproc & nproc are probably not useful)
codes     = names(cancors)[c(9:21)[-c(1,9)]]


############################ Label data: CanCORS and CRN
baseline       = c("ageatdx","stagecat","female")
tmp            = c("case","rel_mths_surg_recurrdt","analysisfu","anyrecurr",
                   "month",codes,baseline)
cancors        = cancors[,tmp]
crn            = crn[,tmp]
CanCORS_PatNum = unique(cancors$case)
CRN_PatNum     = unique(crn$case)
nn             = length(CRN_PatNum)
nnv            = length(CanCORS_PatNum)


CRN_Surv    = data.frame(PatNum = CRN_PatNum, 
                         Delta = sapply(CRN_PatNum,function(x) {
                           tmp = which(crn$case==x)
                           as.numeric(crn$anyrecurr[tmp[1]])
                         }, USE.NAMES = FALSE), 
                         stringsAsFactors = FALSE)
CRN_Surv$SX = with(crn,pmin(analysisfu,rel_mths_surg_recurrdt))[match(CRN_PatNum,crn$case)]
CRN_Surv$SC = crn$analysisfu[match(CRN_PatNum,crn$case)]
## CanCORS have a shorter follow-up, truncate CRN
CRN_Surv$Delta = as.numeric(with(CRN_Surv, Delta & SC<=100))
CRN_Surv$SC = pmin(CRN_Surv$SC,100)
CRN_Surv$SX = with(CRN_Surv,pmin(SX,SC))

tmp2        = crn[!duplicated(crn$case),c("case",baseline)]
CRN_Surv    = cbind(CRN_Surv,tmp2[match(CRN_Surv$PatNum,tmp2$case),-1])
for(i in 2:3){
  txt = paste0("CRN_Surv$StgCat",i,"=as.numeric(CRN_Surv$stagecat==",i,")")
  eval(parse(text=txt))
} 
CRN_Surv    = CRN_Surv[,-6]
crn         = crn[,tmp[!tmp%in%c("anyrecurr","rel_mths_surg_recurrdt",baseline)]] # delete anyrecurr, survival time, and baseline 
crn[crn$month==0,codes] = 0

seer = SEER_med[,tmp[!tmp%in%c("anyrecurr","rel_mths_surg_recurrdt",baseline)]]
seer[seer$month==0,codes] = 0

### CanCORS validation
CanCORS_Surv    = data.frame(PatNum = CanCORS_PatNum, 
                             Delta = sapply(CanCORS_PatNum,function(x) {
                               tmp = which(cancors$case==x)
                               as.numeric(cancors$anyrecurr[tmp[1]])
                             }, USE.NAMES = FALSE), 
                             stringsAsFactors = FALSE)
CanCORS_Surv$SX = with(cancors,pmin(analysisfu,rel_mths_surg_recurrdt))[match(CanCORS_PatNum,cancors$case)]
CanCORS_Surv$SC = cancors$analysisfu[match(CanCORS_PatNum,cancors$case)]

tmp2            = cancors[!duplicated(cancors$case),c("case",baseline)]
CanCORS_Surv    = cbind(CanCORS_Surv,tmp2[match(CanCORS_Surv$PatNum,tmp2$case),-1])
for(i in 2:3){
  txt = paste0("CanCORS_Surv$StgCat",i,"=as.numeric(CanCORS_Surv$stagecat==",i,")")
  eval(parse(text=txt))
} 
CanCORS_Surv    = CanCORS_Surv[,-6]
cancors      = cancors[,tmp[!tmp%in%c("anyrecurr","rel_mths_surg_recurrdt",baseline)]] # delete anyrecurr, survival time, and baseline 
cancors[cancors$month==0,codes] = 0


write.csv(crn,paste0(wkdir2,"crn.csv"),row.names = FALSE)
write.csv(cancors,paste0(wkdir2,"cancors.csv"),row.names = FALSE)
write.csv(seer,paste0(wkdir2,"seer.csv"),row.names = FALSE)
write.csv(CRN_Surv,paste0(wkdir2,"CRN_Surv.csv"),row.names = FALSE)
write.csv(CanCORS_Surv,paste0(wkdir2,"CanCORS_Surv.csv"),row.names = FALSE)


### CRN&CanCORS are too different; 
### mix them and random sample training vs validation

nn  = 600 ## labeled (training) size 
nnv = nrow(CanCORS_Surv)+nrow(CRN_Surv)-nn ## validation size

TrainSurv   = rbind(CanCORS_Surv,CRN_Surv)
set.seed(1234)
ValidSample = sample(1:(nn+nnv),nnv)
ValidSurv   = TrainSurv[sort(ValidSample),]
TrainSurv   = TrainSurv[-ValidSample,]


TrainCode   = rbind(cancors,crn,seer)
ValidSample = TrainCode$case%in%ValidSurv$PatNum
ValidCode   = TrainCode[ValidSample,]
TrainCode   = TrainCode[!ValidSample,]

TrainPatNum = unique(TrainCode$case)
ValidPatNum = unique(ValidCode$case)
TrainN = cbind(aggregate(TrainCode[,codes[1]],by=list(TrainCode$case),FUN=sum),
               sapply(2:length(codes),function(i) aggregate(TrainCode[,codes[i]],
                                                            by=list(TrainCode$case),FUN=sum)[,2]))
TrainN = TrainN[match(TrainPatNum,TrainN[,1]),]
colnames(TrainN) = c("PatientNum",paste0(codes,"_total"))

ValidN           = cbind(aggregate(ValidCode[,codes[1]],by=list(ValidCode$case),FUN=sum),
                         sapply(2:length(codes),function(i) aggregate(ValidCode[,codes[i]],
                                                                      by=list(ValidCode$case),FUN=sum)[,2]))
ValidN           = ValidN[match(ValidPatNum,ValidN[,1]),]
colnames(ValidN) = c("PatientNum",paste0(codes,"_total"))


write.csv(TrainSurv,paste0(wkdir2,"TrainSurv.csv"),row.names = FALSE)
write.csv(ValidSurv,paste0(wkdir2,"ValidSurv.csv"),row.names = FALSE)

write.csv(TrainCode,paste0(wkdir2,"TrainCode.csv"),row.names = FALSE)
write.csv(ValidCode,paste0(wkdir2,"ValidCode.csv"),row.names = FALSE)

write.csv(TrainN,paste0(wkdir2,"TrainN.csv"),row.names = FALSE)
write.csv(ValidN,paste0(wkdir2,"ValidN.csv"),row.names = FALSE)



