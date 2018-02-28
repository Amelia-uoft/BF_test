
args <- commandArgs(TRUE)


region = args[1]
ss = eval( parse(text=args[2]) )
jj = eval( parse(text=args[3]) )

filename = list.files("/home/briollaislab/jxu/BayesFactor/sim1000G/VCF_long")
# file.size = function(file){
#   read.
#   return(file.info(file)$size)
# }
# filename = filename[sapply(filename, file.size) > 25000] #>25kb


if(region=="50k"){
  
  # sample.size = ss*100/2
  nvariants = 72
  p=8
}else if(region=="100k"){
  
  # sample.size = ss*100/2
  nvariants = 147
  p=15
}else if(region=="300k"){
  
  # sample.size = ss*100/2
  nvariants = 442
  p=45
}

library(sim1000G) # I should not use lib.loc
## Download and use full chromosome genetic map
downloadGeneticMap(4)
readGeneticMap(4)
# simulate dataset under H1 using sim1000G
# install.packages("/home/briollaislab/jxu/R/local.package/sim1000G_1.36.2.tar",repos=NULL)
#   library(sim1000G,lib.loc = "/home/briollaislab/jxu/R/x86_64-pc-linux-gnu-library/3.3")

convert.gt = function(xx){
  xx = ifelse(xx==2,1,xx)
  return(xx)
}

data_sim_H0 = function(vcf_sub,seed.num){
  startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) #what does this step do, usually how can I choose number of individuals
  SIM$reset() #what does this step do 
  
  set.seed(seed.num)
  id = generateUnrelatedIndividuals(sample.size)
  # print(id)
  gt = retrieveGenotypes(id)
  gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
  freq = apply(gt,2,sum)/(2*nrow(gt))
  if(sum(freq>=0.5)>0){
    for(jj in which(freq>=0.5)){
      gt[,jj] = 1-gt[,jj]
    }
  }
  if(sum(freq>=0.05)>0){
    gt = gt[,-which(freq>=0.05)]
  }
#  freq = apply(gt,2,sum)/(2*nrow(gt))
  set.seed(seed.num)
  genocase = sample(sample.size,ss)
  return(list(geno = gt[genocase,],nvar = ncol(gt)))
}

data_sim = function(vcf_sub,seed.num){
  startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) #what does this step do, usually how can I choose number of individuals
  SIM$reset() #what does this step do 
  
  set.seed(seed.num)
  id = generateUnrelatedIndividuals(sample.size)
  # print(id)
  gt = retrieveGenotypes(id)
  gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
  freq = apply(gt,2,sum)/(2*nrow(gt))
  if(sum(freq>=0.5)>0){
    for(jj in which(freq>=0.5)){
      gt[,jj] = 1-gt[,jj]
    }
  }
  if(sum(freq>=0.05)>0){
    gt = gt[,-which(freq>=0.05)]
  }
  # #  freq = apply(gt,2,sum)/(2*nrow(gt))
  # set.seed(seed.num)
  # genocase = sample(sample.size,ss)
  # return(list(geno = gt[genocase,],nvar = ncol(gt)))
  freq = apply(gt,2,sum)/(2*nrow(gt))
  causal = sample(setdiff(1:ncol(gt),which(freq==0)),p)
  
  beta.sign = rep(1,p)
  c.value = 0.402
  beta.abs = c.value*abs(log10(freq[causal]))
  beta.val = beta.sign*beta.abs
  x.bar = apply(gt[,causal],2,mean)
  x.bar = as.matrix(x.bar)
  beta.val = t(as.matrix(beta.val))
  beta0 = 0-beta.val %*% x.bar
  
  eta = beta.val %*% t(gt[,causal])
  eta = as.vector(eta) + rep(beta0,nrow(gt))
  prob = exp(eta)/(1+exp(eta))
  
  genocase = rep(NA, sample.size)
  # set.seed(seed.num)
  for(i in 1:sample.size){
    genocase[i] = rbinom(1, 1, prob[i])
  }
  case.idx = sample(which(genocase==1),ss/2)
  control.idx = sample(which(genocase==0),ss/2)
  return(geno = gt[c(case.idx,control.idx),])
  
}
source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/Fun_mix_eta_w_09Jan2018.R")
# source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/Fun_mix_EM_05Dec2017.R")
source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/Fun_mix_EM_15Jan2018_marg_lik.R")
# source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/Fun_reg_MLE_11Nov2017.R")
source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/Fun_reg_MLE_K_03Jan2018_marg_lik.R")
source("/home/briollaislab/jxu/BayesFactor/Simulation/Rcode/subset_vcf_11Jan2018.R")
result.reg= NULL
result.mix= NULL
result.mix.both = NULL


  # for(seed_number in 1:4){ 
  print(filename[jj])
  vcf_file = paste("/home/briollaislab/jxu/BayesFactor/sim1000G/VCF_long/",filename[jj],sep="") 
  vcf = readVCF( vcf_file, maxNumberOfVariants = 1e6 , min_maf = 1e-6 ,max_maf = 0.01)
  
  print(vcf_file)
  for(kk in 1001:1076){
    set.seed(kk)
    
    vcf_sub = subset_vcf(kk)
    sample.size=ss*1.5
    
    # under H0
    simData0 = data_sim_H0(vcf_sub,kk)
    Z0.skat = simData0$geno
    nsites = simData0$nvar
    # Z1.skat = data_sim(vcf_sub,kk)
    # write.csv(Z0.skat,paste("/home/briollaislab/jxu/BayesFactor/sim1000G/result/gtData_VCF",seed_number,"_05Jan2018.csv",sep=""),row.names=F)
    dataset_0 = cbind(apply(Z0.skat,1,sum),c(rep(1,ss/2),rep(0,ss/2)))
    dataset_0 = data.frame(dataset_0)
    # dataset_1 = cbind(apply(Z1.skat,1,sum),c(rep(1,ss/2),rep(0,ss/2)))
    # dataset_1 = data.frame(dataset_1)
    
    result.mix = rbind(result.mix,c(jj,kk,nsites,eta_par_mix(dataset_0,nsites)))
    # result.mix = rbind(result.mix,c(jj,kk,nsites,eta_par_mix(dataset_0,nsites),eta_par_mix(dataset_1,nsites)))
    causal.pool = 1:nsites
    # par.theory0 = eta_par_reg(dataset_0,500)
    # par.theory1 = eta_par_reg(dataset_1,500)
    par.theory0 = eta_par_reg_K(dataset_0)
    # par.theory1 = eta_par_reg_K(dataset_1)
    
    # result.reg = rbind(result.reg,c(jj,kk,nsites,par.theory0,par.theory1))
    result.reg = rbind(result.reg,c(jj,kk,nsites,par.theory0))
    result.mix = data.frame(result.mix)
    result.reg = data.frame(result.reg)
    
    # result.mix.both = rbind(result.mix.both,c(jj,kk,nsites,eta_w_mix(dataset_0,nsites),eta_w_mix(dataset_1,nsites)))
    # result.mix.both = data.frame(result.mix.both)

    # 
    # names(result.mix.both) = c("VCF_file","seed_number","sites",
    #                            "K_0","eta.tilde.total_0", "k.tilde.total_0"," eta.tilde.case_0", "k.tilde.case_0","eta.tilde.control_0", "k.tilde.control_0",
    #                            "eta_total_0","var_total_0","eta_case_0","var_case_0","eta_control_0","var_control_0","BF0",
    #                            "K_1","eta.tilde.total_1", "k.tilde.total_1"," eta.tilde.case_1", "k.tilde.case_1","eta.tilde.control_1", "k.tilde.control_1",
    #                            "eta_total_1","var_total_1","eta_case_1","var_case_1","eta_control_1","var_control_1","BF1")
    
    names(result.mix) = c("VCF_file","seed_number","sites"
                          ,"w0_0","K","eta_total_0","eta_case_0","eta_control_0","loglike_0","loglike1_0","loglike2_0","BF0")

    # names(result.mix) = c("VCF_file","seed_number","sites","w0_0","eta_total_0","var_total_0","eta_case_0","var_case_0","eta_control_0","var_control_0","converge_0","BF0"
    #                       ,"w0_1","eta_total_1","var_total_1","eta_case_1","var_case_1","eta_control_1","var_control_1","converge_1","BF1")
    names(result.reg) = c("VCF_file","seed_number","sites","K_0","eta_total_0",
                          "eta_case_0","eta_control_0","loglike_0","loglike1_0","loglike2_0","BF0")
  }

write.csv(result.mix,paste("/home/briollaislab/jxu/BayesFactor/sim1000G/result/BF_mix_IND_",ss,"_",region,"VCF",jj,"_H0_19Jan2018_9.csv",sep=""),row.names = F)
write.csv(result.reg,paste("/home/briollaislab/jxu/BayesFactor/sim1000G/result/BF_reg_IND_",ss,"_",region,"VCF",jj,"_H0_19Jan2017_9.csv",sep=""),row.names = F)
  # write.csv(result.reg,paste("/home/briollaislab/jxu/BayesFactor/sim1000G/result/BF_reg_IND_",ss,"_",region,"VCF",jj,"_11Jan2017_K.csv",sep=""),row.names = F)
# write.csv(data.frame(result.mix.both),paste("/home/briollaislab/jxu/BayesFactor/sim1000G/result/BF_mix_both_IND_",ss,"_",region,"VCF",jj,"_22Jan2018.csv",sep=""),row.names = F)


