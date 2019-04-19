---
title: "통계 그래픽스 HW#3"
output: word_document
---

```{r, include=FALSE}
setwd("C:/Users/user/misaeng/수업/2019-1/통계그래픽스/R코드&과제")
```

```{r, message=FALSE, warning=FALSE}
library(tidyverse)
data <- read_csv("HW3data.csv")
data
```

##### Missing Values 처리
```{r, message=FALSE, warning=FALSE}
summary(data)
apply(data[,3:16], 2, function(x){ sum(x=="<0.500", na.rm=T) })

data$`002_A2_V2(0wk)` <- as.numeric(data$`002_A2_V2(0wk)`)
data$`002_A2_V4(8wk)` <- as.numeric(data$`002_A2_V4(8wk)`)
data$`002_A2_tread_after_V2(0wk)` <- as.numeric(data$`002_A2_tread_after_V2(0wk)`)
data$`002_A2_tread_after_V4(8wk)` <- as.numeric(data$`002_A2_tread_after_V4(8wk)`)

apply(data[,3:16], 2, function(x){ sum(x=="<0.500", na.rm=T) })
summary(data)
```
"<0.500"의 값이 포함된 변수의 경우 type이 character로 지정되어있기 때문에,
type을 numeric으로 바꿔주면 "<0.500"와 같은 값들은 NA로 처리된다.


##### gather & spread 함수를 이용하여 data 형식 바꾸기
```{r, message=FALSE, warning=FALSE}
data2 <- data %>% gather(3:18, key="Measurement", value="Score") %>%
  mutate(Week = factor(c(rep(rep(c("-2week", "0week", "8week"), each=80), 4), rep(rep(c("0week", "8week"), each=80), 2)))) %>%
  select(-3) %>% mutate(Measurement = factor(c(rep(c("Var_A1", "Var_A1_tread_after", "002_A2", "002_A2_tread_after"), each=240), rep(c("Var_A3", "Var_A3_tread_after"), each=160)), levels=c("Var_A1", "Var_A1_tread_after", "002_A2", "002_A2_tread_after", "Var_A3", "Var_A3_tread_after")))

data2
data2 %>% spread(key=Week, value=Score)
```

##### dplyr과 ggplot을 이용한 EDA
```{r, message=FALSE, warning=FALSE}
Mean <- group_by(data2, Group, Measurement, Week) %>% summarise(Mean_Score = mean(Score, na.rm=T))

ggplot(Mean %>% filter(Measurement=="Var_A1")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("Var_A1 측정 결과 (평균값)")
ggplot(Mean %>% filter(Measurement=="Var_A1_tread_after")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("운동 후 Var_A1 재측정 결과 (평균값)")
ggplot(Mean %>% filter(Measurement=="002_A2")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("002_A2 측정 결과 (평균값)")
ggplot(Mean %>% filter(Measurement=="002_A2_tread_after")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("운동 후 002_A2 재측정 결과 (평균값)")
ggplot(Mean %>% filter(Measurement=="Var_A3")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("Var_A3 측정 결과 (평균값)")
ggplot(Mean %>% filter(Measurement=="Var_A3_tread_after")) + geom_col(aes(x=Group, y=Mean_Score, fill=Group)) +
  facet_grid(. ~ Week) + ggtitle("운동 후 Var_A3 재측정 결과 (평균값)")
```
