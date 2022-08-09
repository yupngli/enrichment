#######################################################################
#
# Pfizer- Sequential Enrichment Designs for Early Phase Clinical Trials
#
#######################################################################

library(dplyr)
library(ggplot2)
library(Rcpp)
library(reshape2)
library(stargazer) # 给txt输出排版用的包

#setwd("C:/Users/yli01/OneDrive - lianbio/working/Enrichment/enrichment/序贯富集_多阶段序贯设计")

start_time <- Sys.time()

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 0: 导入计算后验的c++外部函数
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

sourceCpp('function_cpp.cpp')
sourceCpp('predprobMedian.cpp')
sourceCpp('predprobPosterior.cpp')

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 主程序1: 寻找中期分析的决策边界
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

# 先验
beta_a = 0.1904762 #先验变化会对minimum sample size有影响，就算是无信息先验
beta_b = 1

# 双标准决策中的边界
postmed = 0.24
nv = 0.16
clinicalsig = 0.75
pp_thres = 0.1

# 搜索的responders上限
max_responders <- 20

# 第一阶段full人群的样本量和final阶段full人群样本量,用来算PP
full_n1 <- 50
full_n_total <- 100
# 第一阶段Y+人群的样本量和final阶段Y+人群样本量,用来算PP
y_plus_n1 <- 30
y_plus_n_total <- 60
# 第一阶段Y-人群的样本量和final阶段Y-人群样本量,用来算PP
y_minus_n1 <- 20
y_minus_n_total <- 40

# 设置向量来收集full/Y+/Y-在期中的时候PP的结果
r_vec <- seq(1,max_responders,1)
full_succ_vec <- rep(NA,max_responders)
y_minus_succ_vec <- rep(NA,max_responders)
y_plus_succ_vec <- rep(NA,max_responders)

for (r in 0:max_responders){
  full_pp <- predprobMedian(r,full_n1,full_n_total,beta_a,beta_b,postmed)
  y_minus_pp <- predprobPosterior(r,y_minus_n1,y_minus_n_total,beta_a,beta_b,nv,clinicalsig)
  y_plus_pp <- predprobMedian(r,y_plus_n1,y_plus_n_total,beta_a,beta_b,postmed)
  
  # collect result
  full_succ_vec[r] <- full_pp
  y_minus_succ_vec[r] <- y_minus_pp
  y_plus_succ_vec[r] <- y_plus_pp
  
  # package result
  result_df <- data.frame(r=r_vec,
                          fullpp=full_succ_vec,
                          yminuspp=y_minus_succ_vec,
                          ypluspp=y_plus_succ_vec)
}

result_df_t <- melt(result_df,id.vars = "r",variable.name = "population_pp")
result_df_t <- result_df_t %>% filter(!is.na(value)) 

p <- ggplot(data=result_df_t, mapping = aes(x=r,y=value)) +
  geom_point(aes(color = factor(population_pp), shape=factor(population_pp)), size = 3) + 
  geom_hline(yintercept=clinicalsig, linetype="dotted", color = "orange",lwd=1) +
  geom_hline(yintercept=pp_thres, linetype="dotted", color = "orange",lwd=1) +
  # geom_vline(xintercept=58, linetype="dashed", color = "black") +
  # geom_vline(xintercept=87, linetype="dashed", color = "black") +
  scale_y_continuous(name="Predictive Probability")+
  scale_x_continuous(name="#responders",breaks = seq(1,max_responders,1))+
  labs(color="Group",shape="Group") + 
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major.x = element_line(colour = "grey70",linetype = "longdash"))
p

# 打印决策边界
result_df_t_decision <- result_df_t %>%
  mutate(gt10=ifelse(value>0.1,TRUE,FALSE),
         gt75=ifelse(value>0.1,TRUE,FALSE))

### Continue with F:POS(F) ≥ 10% and POS(Y-) ≥ 10%,能同时满足这两个条件
### 从上面的图中可以看出，当full人群的中期条件满足（50人里观察到10个full人群responders）,此时去检查y-人群，responders>=3时可以满足Full期中的决策边界
F_dec <- result_df_t_decision %>%
  filter(population_pp %in% c("fullpp","yminuspp") & value >= 0.1) %>%
  group_by(population_pp) %>% 
  slice(which.min(r))
F_dec_r <- max(F_dec$r)
cat("The Go decision for full population is responders >= ",F_dec_r)

### Continue with Y+
### 从上面的图中可以看出，当full人群的中期条件满足（50人里观察到10个full人群responders）,此时去检查y-人群，responders>=3时可以满足Full期中的决策边界
y_plus_dec <- result_df %>% 
  mutate(cond1 = ifelse(fullpp>=0.1 & yminuspp<0.1 & ypluspp>=0.1,TRUE,FALSE),
         cond2 = ifelse(fullpp<0.1 & ypluspp>=0.1,TRUE,FALSE),
         overall = cond1 | cond2)

y_plus_dec2 <- y_plus_dec %>%
  filter(overall == TRUE) %>%
  slice(which.min(r))

y_plus_dec_r <- max(y_plus_dec2$r)
cat("The Go decision for Y+ population is responders >= ",y_plus_dec_r)


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# 主程序2: 评估设计的OC
#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

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

# 人数假设
allcomers_ia <- full_n1
y_plus_ia <- y_plus_n1
y_minus_ia <- y_minus_n1
allcomers_overall <- full_n_total
y_plus_overall <- y_plus_n_total
y_minus_overall <- y_minus_n_total
rep <- 5000

simseqoc <- function(orr_plus,
                  orr_minus,
                  allcomers_ia,
                  y_plus_ia,
                  y_minus_ia,
                  allcomers_overall,
                  y_plus_overall,
                  y_minus_overall,
                  rep){
  
  # 设置responder threshold,从上面的试验设计中得来
  full_IA_gate <- 10
  full_FA_gate <- 25
  
  y_plus_IA_gate <- 6
  y_plus_FA_gate <- 15
  
  y_minus_IA_gate <- 3
  y_minus_FA_gate <- 9
  
  # 初始化模拟矩阵, s1: stage 1; s2: stage 2; overall=s1+s2
  simmat <- data.frame(orr_plus=rep(NA,rep),
                       orr_minus=rep(NA,rep),
                       allcomers_s1=rep(NA,rep),
                       y_plus_s1=rep(NA,rep),
                       y_minus_s1=rep(NA,rep),
                       allcomers_s2=rep(NA,rep),
                       y_plus_s2=rep(NA,rep),
                       y_minus_s2=rep(NA,rep),
                       full_gate_s1=rep(NA,rep),
                       y_plus_gate_s1=rep(NA,rep),
                       y_minus_gate_s1=rep(NA,rep),
                       full_gate_s2=rep(NA,rep),
                       y_plus_gate_s2=rep(NA,rep),
                       y_minus_gate_s2=rep(NA,rep),
                       y_plus_responders_s1=rep(NA,rep),
                       full_responders_s1=rep(NA,rep),
                       y_plus_responders_overall=rep(NA,rep),
                       full_responders_overall=rep(NA,rep),
                       succ_ia_full=rep(FALSE,rep),
                       succ_ia_y_pos=rep(FALSE,rep),
                       early_stop=rep(FALSE,rep),
                       succ_overall_full=rep(FALSE,rep),
                       succ_overall_y_pos=rep(FALSE,rep))
  
  for (row in 1:rep){
    
    # 产生overall样本量的随机数
    y_plus_data <- rbinom(y_plus_overall,1,orr_plus)
    y_minus_data <- rbinom(y_minus_overall,1,orr_minus)
    
    # 切割向量，生成中期分析数据集
    y_plus_ia_data <- y_plus_data[1:y_plus_ia]
    y_minus_ia_data <- y_minus_data[1:y_minus_ia]
    
    # 计算IA 及 overall分别有多少responders
    ## IA
    resp_y_plus_ia <- sum(y_plus_ia_data)
    resp_y_minus_ia <- sum(y_minus_ia_data)
    total_resp_ia <- resp_y_plus_ia+resp_y_minus_ia
    ## OVERALL
    resp_y_plus_overall <- sum(y_plus_data)
    resp_y_minus_overall <- sum(y_minus_data)
    total_resp_overall <- resp_y_plus_overall+resp_y_minus_overall
    
    #决策 - 怎么避免重复书写? if 太多了
    ## 先设置不同IA前进方向/最终结果判断的flag,否则IF太多了很乱
    IA_go_withF <- FALSE
    IA_go_withYPlus <- FALSE
    IA_ET <- FALSE
    SUCC_F <- FALSE
    SUCC_YPlus <- FALSE
    FA_FAIL <- FALSE
    ## S1的决策部分
    if (total_resp_ia>=10){
      if (resp_y_minus_ia>=3){
        IA_go_withF <- TRUE
      } else if (resp_y_minus_ia<3){
        if (resp_y_plus_ia>=6){
          IA_go_withYPlus <- TRUE
        } else {
          IA_ET <- TRUE
        }
      }
    } else if (total_resp_ia<10){
      if (resp_y_plus_ia>=6){
        IA_go_withYPlus <- TRUE
      } else {
        IA_ET <- TRUE
      }
    }
    ## S2的决策部分
    ### Full
    if (IA_go_withF==TRUE){
      if (total_resp_overall >= full_FA_gate){
        if (resp_y_minus_overall >= y_minus_FA_gate){
          SUCC_F <- TRUE
        } else {
          if (resp_y_plus_overall >= y_plus_FA_gate){
            SUCC_YPlus <- TRUE
          }
        }
      } else{
        if (resp_y_plus_overall >= y_plus_FA_gate){
          SUCC_YPlus <- TRUE
        }
      }
    }
    ### Y+
    if (IA_go_withF!=TRUE & IA_go_withYPlus==TRUE){
      if (resp_y_plus_overall >= y_plus_FA_gate){
        SUCC_YPlus <- TRUE
      } else {
        FA_FAIL <- TRUE
      }
    }
    
    
    # simmat赋值
    simmat[row,"orr_plus"] <- orr_plus
    simmat[row,"orr_minus"] <- orr_minus
    simmat[row,"allcomers_s1"] <- allcomers_ia
    simmat[row,"y_plus_s1"] <- y_plus_ia
    simmat[row,"y_minus_s1"] <- y_minus_ia
    simmat[row,"allcomers_s2"] <- allcomers_overall
    simmat[row,"y_plus_s2"] <- y_plus_overall
    simmat[row,"y_minus_s2"] <- y_minus_ia
    simmat[row,"full_gate_s1"] <- full_IA_gate
    simmat[row,"y_plus_gate_s1"] <- y_plus_IA_gate
    simmat[row,"y_minus_gate_s1"] <- y_minus_IA_gate
    simmat[row,"full_gate_s2"] <- full_FA_gate
    simmat[row,"y_plus_gate_s2"] <- y_plus_FA_gate
    simmat[row,"y_minus_gate_s2"] <- y_minus_FA_gate
    simmat[row,"y_plus_responders_s1"] <- resp_y_plus_ia
    simmat[row,"full_responders_s1"] <- total_resp_ia
    simmat[row,"y_plus_responders_overall"] <- resp_y_plus_overall
    simmat[row,"full_responders_overall"] <- total_resp_overall
    simmat[row,"succ_ia_full"] <- IA_go_withF
    simmat[row,"succ_ia_y_pos"] <- IA_go_withYPlus
    simmat[row,"early_stop"] <- IA_ET
    simmat[row,"succ_overall_full"] <- SUCC_F
    simmat[row,"succ_overall_y_pos"] <- SUCC_YPlus
      
  }
  
  return(simmat)
}

rep = 5000

scen1 <- simseqoc(orr_plus_1,
               orr_minus_1,
               allcomers_ia,
               y_plus_ia,
               y_minus_ia,
               allcomers_overall,
               y_plus_overall,
               y_minus_overall,
               rep)

scen2 <- simseqoc(orr_plus_2,
               orr_minus_2,
               allcomers_ia,
               y_plus_ia,
               y_minus_ia,
               allcomers_overall,
               y_plus_overall,
               y_minus_overall,
               rep)

scen3 <- simseqoc(orr_plus_3,
               orr_minus_3,
               allcomers_ia,
               y_plus_ia,
               y_minus_ia,
               allcomers_overall,
               y_plus_overall,
               y_minus_overall,
               rep)

scen4 <- simseqoc(orr_plus_4,
               orr_minus_4,
               allcomers_ia,
               y_plus_ia,
               y_minus_ia,
               allcomers_overall,
               y_plus_overall,
               y_minus_overall,
               rep)

# 汇报结果 - 与壁报的Figure 6.2 OC for One stage Enrichment Design 对比
result <- data.frame(scenario=rep(NA,4),
                     prob_eff_full_IA=rep(NA,4),
                     prob_eff_positive_IA=rep(NA,4),
                     prob_ET_IA=rep(NA,4),
                     prob_eff_full_FA=rep(NA,4),
                     prob_eff_positive_FA=rep(NA,4))

result[1,"scenario"] <- "ORR (Y+) = ORR (Y-) = 16%"
result[2,"scenario"] <- "ORR (Y+) = 32%, ORR (Y-) = 16%"
result[3,"scenario"] <- "ORR (Y+) = 16%, ORR (Y-) = 32%"
result[4,"scenario"] <- "ORR (Y+) = ORR (Y-) = 32%"


result[1,"prob_eff_full_IA"] <- mean(scen1$succ_ia_full)*100
result[2,"prob_eff_full_IA"] <- mean(scen2$succ_ia_full)*100
result[3,"prob_eff_full_IA"] <- mean(scen3$succ_ia_full)*100
result[4,"prob_eff_full_IA"] <- mean(scen4$succ_ia_full)*100

result[1,"prob_eff_positive_IA"] <- mean(scen1$succ_ia_y_pos)*100
result[2,"prob_eff_positive_IA"] <- mean(scen2$succ_ia_y_pos)*100
result[3,"prob_eff_positive_IA"] <- mean(scen3$succ_ia_y_pos)*100
result[4,"prob_eff_positive_IA"] <- mean(scen4$succ_ia_y_pos)*100

result[1,"prob_ET_IA"] <- mean(scen1$early_stop)*100
result[2,"prob_ET_IA"] <- mean(scen2$early_stop)*100
result[3,"prob_ET_IA"] <- mean(scen3$early_stop)*100
result[4,"prob_ET_IA"] <- mean(scen4$early_stop)*100

result[1,"prob_eff_full_FA"] <- mean(scen1$succ_overall_full)*100
result[2,"prob_eff_full_FA"] <- mean(scen2$succ_overall_full)*100
result[3,"prob_eff_full_FA"] <- mean(scen3$succ_overall_full)*100
result[4,"prob_eff_full_FA"] <- mean(scen4$succ_overall_full)*100

result[1,"prob_eff_positive_FA"] <- mean(scen1$succ_overall_y_pos)*100
result[2,"prob_eff_positive_FA"] <- mean(scen2$succ_overall_y_pos)*100
result[3,"prob_eff_positive_FA"] <- mean(scen3$succ_overall_y_pos)*100
result[4,"prob_eff_positive_FA"] <- mean(scen4$succ_overall_y_pos)*100

# 评估程序运行时间
end_time <- Sys.time()

# 打印到txt文件
stargazer(result,
          title = "Interim & Final OC for Sequential Enrichment Design",
          summary = FALSE,
          rownames = FALSE,
          type = "text",
          align = TRUE,
          out = "两阶段富集设计OC.txt",
          notes = paste0("程序运行用时",as.character.Date(end_time-start_time),", 可对比辉瑞壁报Figure 8.2与Figure 8.3"),
          notes.align = "l")