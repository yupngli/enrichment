#######################################################################
#
# Pfizer- One Stage Enrichment Designs for Early Phase Clinical Trials
#
#######################################################################

library(dplyr)
library(ggplot2)
library(Rcpp)
library(stargazer) # 给txt输出排版用的包

start_time <- Sys.time()

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 0: 导入计算后验的c++外部函数
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#setwd("C:/Users/yli01/OneDrive - lianbio/working/Enrichment/enrichment")
sourceCpp('function_cpp.cpp')

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 主程序1: 寻找Minimum sample size - Figure 1
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

nv = 0.16
statsig = 0.975
clinicalsig = 0.75
postmed = 0.24
positivestatsig = 0.95
maxsamplesize = 110 #我们能接受的样本量的上限，但是设太低可能会找不到
samplesize_start = 30 #实际应用中，总不可能样本量从1开始的
maxresponse = 25 #实际应用中可能出现的最高的responders数量?

beta_a = 0.1904762 #先验变化会对minimum sample size有影响，就算是无信息先验
beta_b = 1

totalrow = (maxsamplesize-samplesize_start+1)*(maxresponse)

# bothfl: for Bayesian dual criteria for Y+ patients
# tripfl: for Bayesian triplet criterion for all patients
# df <- data.frame(r=rep(NA,totalrow),
#                  n=rep(NA,totalrow),
#                  significance=rep(NA,totalrow),
#                  postmedian=rep(NA,totalrow),
#                  justmedfl=rep(FALSE,totalrow),
#                  bothfl=rep(FALSE,totalrow),
#                  tripfl=rep(FALSE,totalrow),
#                  group=rep(NA,totalrow))

# row <- 1

# for (r in 1:maxresponse){
#   for (n in samplesize_start:maxsamplesize){
#     if (r>n) {
#       next
#     }
#     rand <- rbeta(100000,beta_a+r,beta_b+n-r)
#     sig <- mean(rand>=nv)
#     med <- median(rand)
#     df[row,1]=r
#     df[row,2]=n
#     df[row,3]=sig
#     df[row,4]=med
#     # 仅满足一个标准(median>xxx)
#     if (sig<positivestatsig & med>=postmed){df[row,5]=TRUE}
#     # 仅满足两个标准(median>xxx与阳性显著)
#     if (sig>=positivestatsig & med>=postmed & sig<statsig){df[row,6]=TRUE}
#     # 满足三个标准(median>xxx与全人群显著)
#     if (sig>=statsig & med>=postmed){df[row,7]=TRUE}
#     # 为绘图分组
#     if (df[row,7]==TRUE){
#       df[row,8] <- "full"
#     } else if(df[row,6]==TRUE){
#       df[row,8] <- "positive only"
#     } else if (df[row,5]==TRUE){
#       df[row,8] <- "N/A"
#     }
#     row <- row+1
#   }
# }

df <- as.data.frame(responsegate(maxresponse = maxresponse,
                   samplesizestart = samplesize_start,
                   maxsamplesize = maxsamplesize,
                   a = beta_a,
                   b = beta_b,
                   nv = nv,
                   positivestatsig = positivestatsig,
                   statsig = statsig,
                   postmed = postmed))

colnames(df) <- c("r","n","significance","postmedian","justmedfl","bothfl","tripfl","group")

# 删去不满足任何标准的行记录
df_cc <- filter(df,justmedfl==1 | bothfl==1 | tripfl==1)

# recode
df_cc <- df_cc %>%
  mutate(group=recode(group,'2'="full",'3'="positive only",'4'="n/a"))

## 最小样本量逻辑：按照n分组，找到最小的那个significance，因为一旦满足当前的要求，
## 则对于这个n，比当前significance对应的r更大的r，它的significance肯定更高

df_cc_plot <- df_cc %>% 
  group_by(n) %>% 
  slice(which.min(significance))

# 可视化最小样本量
text_y <- min(df_cc_plot$significance)+0.002
min_y <- min(df_cc_plot$significance)
max_y <- max(df_cc_plot$significance)

p <- ggplot(data=df_cc_plot, mapping = aes(x=n,y=significance)) +
  geom_point(aes(color = factor(group)), size = 2, shape=17) + 
  geom_hline(yintercept=positivestatsig, linetype="dashed", color = "orange") +
  geom_hline(yintercept=statsig, linetype="dashed", color = "orange") +
  geom_vline(xintercept=58, linetype="dashed", color = "black") +
  geom_vline(xintercept=87, linetype="dashed", color = "black") +
  annotate(geom="text", x=58+7, y=text_y, label= "r/n=15/58") + 
  annotate(geom="text", x=87+7, y=text_y, label= "r/n=22/87") + 
  scale_y_continuous(name="Statistical Significance", limits=c(min_y, max_y))+
  labs(color="Group") + 
  theme_classic()
p

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 主程序2: 评估设计的OC
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

allcomers = 100
y_plus = 60
y_minus = allcomers-y_plus

# 基于100人总样本量，60人Y+，40人Y-,去数据集df_cc_plot里寻找对应的responder gate
responder_gate <- df_cc_plot %>% 
  filter(n == 100 | n==60 | n==40) %>%
  select(c(r,n))

full_gate <- responder_gate[which(responder_gate$n==allcomers),1]$r
y_plus_gate <- responder_gate[which(responder_gate$n==y_plus),1]$r
y_minus_gate <- responder_gate[which(responder_gate$n==y_minus),1]$r

# 不同场景的orr假设
## 1
orr_plus_1 <- 0.16
orr_minus_1 <- 0.16
## 2
orr_plus_2 <- 0.32
orr_minus_2 <- 0.16
## 3
orr_plus_3 <- 0.16
orr_minus_3 <- 0.32
## 4
orr_plus_4 <- 0.32
orr_minus_4 <- 0.32

simoc <- function(orr_plus,
                  orr_minus,
                  allcomers,
                  y_plus,
                  y_minus,
                  full_gate,
                  y_plus_gate,
                  y_minus_gate,
                  rep){
  
  # 初始化模拟矩阵
  simmat <- data.frame(orr_plus=rep(NA,rep),
                   orr_minus=rep(NA,rep),
                   allcomers=rep(NA,rep),
                   y_plus=rep(NA,rep),
                   y_minus=rep(NA,rep),
                   full_gate=rep(NA,rep),
                   y_plus_gate=rep(NA,rep),
                   y_plus_responders=rep(NA,rep),
                   full_responders=rep(NA,rep),
                   succ_full=rep(FALSE,rep),
                   succ_y_pos=rep(FALSE,rep))
  
  for (row in 1:rep){
    resp_y_plus <- sum(rbinom(y_plus,1,orr_plus))
    resp_y_minus <- sum(rbinom(y_minus,1,orr_minus))
    total_resp <- resp_y_minus+resp_y_plus
    
    #决策 - 怎么避免重复书写? if 太多了
    if (total_resp >= full_gate){
      if (resp_y_minus >= y_minus_gate){
        simmat[row,"succ_full"] <- TRUE
      } else {
        if (resp_y_plus >= y_plus_gate){
          simmat[row,"succ_y_pos"] <- TRUE
        }
      }
    } else{
      if (resp_y_plus >= y_plus_gate){
        simmat[row,"succ_y_pos"] <- TRUE
      }
    }
    
    # simmat赋值
    simmat[row,"orr_plus"] <- orr_plus
      simmat[row,"orr_minus"] <- orr_minus
      simmat[row,"allcomers"] <- allcomers
      simmat[row,"y_plus"] <- y_plus
      simmat[row,"y_minus"] <- y_minus
      simmat[row,"full_gate"] <- full_gate
      simmat[row,"y_plus_gate"] <- y_plus_gate
      simmat[row,"y_plus_responders"] <- resp_y_plus
      simmat[row,"full_responders"] <- total_resp
  }
  
  return(simmat)
}

rep = 12000

scen1 <- simoc(orr_plus_1,
               orr_minus_1,
               allcomers,
               y_plus,
               y_minus,
               full_gate,
               y_plus_gate,
               y_minus_gate,
               rep)

scen2 <- simoc(orr_plus_2,
               orr_minus_2,
               allcomers,
               y_plus,
               y_minus,
               full_gate,
               y_plus_gate,
               y_minus_gate,
               rep)

scen3 <- simoc(orr_plus_3,
               orr_minus_3,
               allcomers,
               y_plus,
               y_minus,
               full_gate,
               y_plus_gate,
               y_minus_gate,
               rep)

scen4 <- simoc(orr_plus_4,
               orr_minus_4,
               allcomers,
               y_plus,
               y_minus,
               full_gate,
               y_plus_gate,
               y_minus_gate,
               rep)

# 汇报结果 - 与壁报的Figure 6.2 OC for One stage Enrichment Design 对比
result <- data.frame(scenario=rep(NA,4),
                     prob_eff_full=rep(NA,4),
                     prob_eff_positive=rep(NA,4))

result[1,"scenario"] <- "ORR (Y+) = ORR (Y-) = 16%"
result[2,"scenario"] <- "ORR (Y+) = 32%, ORR (Y-) = 16%"
result[3,"scenario"] <- "ORR (Y+) = 16%, ORR (Y-) = 32%"
result[4,"scenario"] <- "ORR (Y+) = ORR (Y-) = 32%"


result[1,"prob_eff_full"] <- mean(scen1$succ_full)*100
result[2,"prob_eff_full"] <- mean(scen2$succ_full)*100
result[3,"prob_eff_full"] <- mean(scen3$succ_full)*100
result[4,"prob_eff_full"] <- mean(scen4$succ_full)*100

result[1,"prob_eff_positive"] <- mean(scen1$succ_y_pos)*100
result[2,"prob_eff_positive"] <- mean(scen2$succ_y_pos)*100
result[3,"prob_eff_positive"] <- mean(scen3$succ_y_pos)*100
result[4,"prob_eff_positive"] <- mean(scen4$succ_y_pos)*100

# 评估程序运行时间
end_time <- Sys.time()
end_time-start_time

# 输出TXT
stargazer(result,
          title = "OC for One stage Enrichment Design",
          summary = FALSE,
          rownames = FALSE,
          type = "text",
          align = TRUE,
          out = "单阶段富集设计OC.txt",
          notes = paste0("程序运行用时",as.character.Date(end_time-start_time),", 可对比辉瑞壁报Figure 6.2"),
          notes.align = "l")