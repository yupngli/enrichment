library(Rcpp)
sourceCpp('function_cpp.cpp')

# Rcpp

start_time <- Sys.time()
aaa <- responsegate(maxresponse = 100,samplesizestart = 30,maxsamplesize = 200,a = 0.1904762,b = 1,nv = 0.16,positivestatsig = 0.95,statsig = 0.975,postmed = 0.24)
end_time <- Sys.time()

# R
start_time2 <- Sys.time()
nv = 0.16
statsig = 0.975
clinicalsig = 0.75
postmed = 0.24
positivestatsig = 0.95
maxsamplesize = 200 #我们能接受的样本量的上限，但是设太低可能会找不到
samplesize_start = 30 #实际应用中，总不可能样本量从1开始的
maxresponse = 100 #实际应用中可能出现的最高的responders数量?

beta_a = 0.1904762 #先验变化会对minimum sample size有影响，就算是无信息先验
beta_b = 1

totalrow = (maxsamplesize-samplesize_start+1)*(maxresponse)

# bothfl: for Bayesian dual criteria for Y+ patients
# tripfl: for Bayesian triplet criterion for all patients
df <- data.frame(r=rep(NA,totalrow),
                 n=rep(NA,totalrow),
                 significance=rep(NA,totalrow),
                 postmedian=rep(NA,totalrow),
                 justmedfl=rep(FALSE,totalrow),
                 bothfl=rep(FALSE,totalrow),
                 tripfl=rep(FALSE,totalrow),
                 group=rep(NA,totalrow))

row <- 1

for (r in 1:maxresponse){
  for (n in samplesize_start:maxsamplesize){
    if (r>n) {
      next
    }
    rand <- rbeta(100000,beta_a+r,beta_b+n-r)
    sig <- mean(rand>=nv)
    med <- median(rand)
    df[row,1]=r
    df[row,2]=n
    df[row,3]=sig
    df[row,4]=med
    # 仅满足一个标准(median>xxx)
    if (sig<positivestatsig & med>=postmed){df[row,5]=TRUE}
    # 仅满足两个标准(median>xxx与阳性显著)
    if (sig>=positivestatsig & med>=postmed & sig<statsig){df[row,6]=TRUE}
    # 满足三个标准(median>xxx与全人群显著)
    if (sig>=statsig & med>=postmed){df[row,7]=TRUE}
    # 为绘图分组
    if (df[row,7]==TRUE){
      df[row,8] <- "full"
    } else if(df[row,6]==TRUE){
      df[row,8] <- "positive only"
    } else if (df[row,5]==TRUE){
      df[row,8] <- "N/A"
    }
    row <- row+1
  }
}
end_time2 <- Sys.time()

# compare
cat("Rcpp takes: ",end_time - start_time," secs")
cat("R code takes: ",end_time2 - start_time2," secs")
