load("C:/Users/Jayaditya Nath/Documents/old.machine.rda")
data_old = data.frame(old.machine)
data_new = data.frame(new.machine)
data_old
#Plotting the entire data
pairs(data_old)

#Using OLS to approximate the starting values, 

mod_ols = lm(log(Y) ~ log(X1) + log(X2), data = data_old)
coef_ols = coef(mod_ols)


#Initial theta values
initial_theta_curl = c(exp(coef_ols[1]),coef_ols[2],coef_ols[3])
initial_theta_curl


#Given model
model = function(theta_vec,x_1,x_2)
{
  theta_0 = theta_vec[1]
  theta_1 = theta_vec[2]
  theta_2 = theta_vec[3]
  y = theta_0*((x_1)^theta_1)*((x_2)^theta_2)
  return(y)
}



#Gauss_Newton Algorithm to find the final parameter estimates

library(pracma)
Y = as.matrix(data_old$Y,nrow=18)
theta_curl = initial_theta_curl

for(iter in 1:20)
{
  f_theta = as.matrix(model(theta_curl,data_old$X1,data_old$X2),nrow=18)
  col_1 = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_0'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x_1=data_old$X1,x_2=data_old$X2))
  col_2 = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_1'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x_1=data_old$X1,x_2=data_old$X2))
  col_3 = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_2'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x_1=data_old$X1,x_2=data_old$X2))
  F_theta = cbind(col_1,col_2,col_3)
  
  delta = pinv(t(F_theta) %*% F_theta) %*% t(F_theta) %*% (Y - f_theta)
  theta_curl = as.matrix(theta_curl) + delta
  
  if(max(abs(delta))<0.0001)
  {
    break
  }
}


print(paste0("The final estimates of theta are : ",theta_curl[1],",",theta_curl[2]," and ",theta_curl[3]," respectively."))
print(paste0("The algorithm converges at the ",iter," th iteration."))


#Using nls() to check the final estimates found after using the Gauss_Newton Algorithm
mod_nls = nls(Y ~ theta0 * X1^theta1 * X2^theta2, data = data_old, start = list(theta0 = initial_theta_curl[1], theta1 = initial_theta_curl[2], theta2 = initial_theta_curl[3]))
summary(mod_nls)

install.packages("marginaleffects")
library(marginaleffects)
library(car)
deltaMethod(mod_nls,"theta1/theta2")
hypotheses(mod_nls,"theta0/theta1=1")

r_squared_pseudo_old = 1 - (sum(((data_old$Y) - (as.matrix(model(theta_curl,(data_old$X1),(data_old$X2)),nrow=18)))^2)/sum(((data_old$Y) - mean(data_old$Y))^2))

sigma = sqrt((t(Y - f_theta) %*% (Y - f_theta))/(nrow(data_old)-length(theta_curl)))
print(paste0("95% confidence interval for theta_0 is (",theta_curl[1]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),",",theta_curl[1]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),")"))
print(paste0("95% confidence interval for theta_0 is (",theta_curl[2]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),",",theta_curl[2]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),")"))
print(paste0("95% confidence interval for theta_0 is (",theta_curl[3]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),",",theta_curl[3]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta) %*% F_theta)[2,2])),")"))

#Considering the residuals
install.packages("nlstools")
library(nlstools)


par(mfrow=c(2,2))
plot(predict(mod_nls),residuals(mod_nls),xlab = "Fitted values",ylab = "Residuals",main = "Old Machine")
abline(0,0,col="red")
plot(predict(mod_nls_new),residuals(mod_nls_new),xlab = "Fitted values",ylab = "Residuals",main = "New Machine")
abline(0,0,col="red")
qqnorm(data_old$Y-predict(mod_nls))
qqline(data_old$Y-predict(mod_nls),col="red")
qqnorm(data_new$Y-predict(mod_nls_new))
qqline(data_new$Y-predict(mod_nls_new),col="red")

#For the new machine 


#Plotting the entire data
pairs(data_new)

#Using OLS to approximate the starting values, 

mod_ols_new = lm(log(Y) ~ log(X1) + log(X2), data = data_new)
coef_ols_new = coef(mod_ols_new)


#Initial theta values
initial_theta_curl_new = c(exp(coef_ols_new[1]),coef_ols_new[2],coef_ols_new[3])
initial_theta_curl_new


#Gauss_Newton Algorithm to find the final parameter estimates

library(pracma)
Y_new = as.matrix(data_new$Y,nrow=18)
theta_curl_new = initial_theta_curl_new

for(iter in 1:10)
{
  f_theta_new = as.matrix(model(theta_curl_new,data_new$X1,data_new$X2),nrow=18)
  col_1_new = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_0'),list(theta_0=theta_curl_new[1],theta_1=theta_curl_new[2],theta_2=theta_curl_new[3],x_1=data_new$X1,x_2=data_new$X2))
  col_2_new = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_1'),list(theta_0=theta_curl_new[1],theta_1=theta_curl_new[2],theta_2=theta_curl_new[3],x_1=data_new$X1,x_2=data_new$X2))
  col_3_new = eval(D(expression(theta_0*((x_1)^theta_1)*((x_2)^theta_2)),'theta_2'),list(theta_0=theta_curl_new[1],theta_1=theta_curl_new[2],theta_2=theta_curl_new[3],x_1=data_new$X1,x_2=data_new$X2))
  F_theta_new = cbind(col_1_new,col_2_new,col_3_new)
  
  delta_new = pinv(t(F_theta_new) %*% F_theta_new) %*% t(F_theta_new) %*% (Y_new - f_theta_new)
  theta_curl_new = as.matrix(theta_curl_new) + delta_new
  
  if(max(abs(delta_new))<0.0001)
  {
    break
  }
}


print(paste0("The final estimates of theta are : ",theta_curl_new[1],",",theta_curl_new[2]," and ",theta_curl_new[3]," respectively."))
print(paste0("The algorithm converges at the ",iter," th iteration."))


#Using nls() to check the final estimates found after using the Gauss_Newton Algorithm
mod_nls_new =  nls(Y ~ theta0 * X1^theta1 * X2^theta2, data = data_new, start = list(theta0 = initial_theta_curl_new[1], theta1 = initial_theta_curl_new[2], theta2 = initial_theta_curl_new[3]))
summary(mod_nls_new)


#Considering the residuals
plot(nlsResiduals(mod_nls_new))


#Comparing the two models 
anova(mod_nls,mod_nls_new)
predict(mod_nls)
predict(mod_nls_new)
plot(data_old$Y, predict(mod_nls), col = "blue", pch = 16, xlab = "Observed", ylab = "Predicted", main = "Model Comparison")
points(data_new$Y, predict(mod_nls_new), col = "red", pch = 16)
legend("topright", legend = c("Model Old", "Model New"), col = c("blue", "red"), pch = 16)
install.packages("lmtest")
library(lmtest)
lrtest(mod_nls,mod_nls_new)
qqnorm(data_old$Y-predict(mod_nls))
qqline(data_old$Y-predict(mod_nls))
res_old = data_old$Y-predict(mod_nls)
res_new = data_new$Y-predict(mod_nls_new)
data.frame(Residuals_Old = res_old)
cbind(data.frame(Residuals_New=res_new),data.frame(Residuals_Old = res_old))


 #Finding confidence interval for parameter estimates 

sigma = sqrt((t(Y_new - f_theta_new) %*% (Y_new - f_theta_new))/(nrow(data_new)-length(theta_curl_new)))
print(paste0("95% confidence interval for theta_0 is (",theta_curl_new[1]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),",",theta_curl_new[1]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),")"))
print(paste0("95% confidence interval for theta_0 is (",theta_curl_new[2]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),",",theta_curl_new[2]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),")"))
print(paste0("95% confidence interval for theta_0 is (",theta_curl_new[3]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),",",theta_curl_new[3]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_new) %*% F_theta_new)[2,2])),")"))




























library(gridExtra)

set.seed(153535)
data_2 = data.frame(read.csv("C:/Users/Jayaditya Nath/Documents/Q2.csv",header = T))
View(data_2)
dim(data_2)
data.training.id <- sample(1:dim(data_2)[1], 4500, replace = FALSE)
data_training <- data_2[data.training.id,]
data_test <- data_2[-data.training.id,]
dim(data_training)
dim(data_test)
colnames(data_2)
#EDA

pairs(data_training)
cor(data_training)
#No issue of multicollinearity

library(gridExtra)
library(ggplot2)

#Age with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$Age, data_training$Gender, data_training$CHF)
interaction.plot(data_training$Age, data_training$CAD, data_training$CHF)
interaction.plot(data_training$Age, data_training$Diabetes, data_training$CHF)




#BMI with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$BMI, data_training$Gender, data_training$CHF)
interaction.plot(data_training$BMI, data_training$CAD, data_training$CHF)
interaction.plot(data_training$BMI, data_training$Diabetes, data_training$CHF)


#Race with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$Race, data_training$Gender, data_training$CHF)
interaction.plot(data_training$Race, data_training$CAD, data_training$CHF)
interaction.plot(data_training$Race, data_training$Diabetes, data_training$CHF)


#Chol with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$Chol, data_training$Gender, data_training$CHF)
interaction.plot(data_training$Chol, data_training$CAD, data_training$CHF)
interaction.plot(data_training$Chol, data_training$Diabetes, data_training$CHF)


#PA with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$PA, data_training$Gender, data_training$CHF)
interaction.plot(data_training$PA, data_training$CAD, data_training$CHF)
interaction.plot(data_training$PA, data_training$Diabetes, data_training$CHF)


#Alcohol with CHF
par(mfrow=c(1,3))
interaction.plot(data_training$Alcohol, data_training$Gender, data_training$CHF)
interaction.plot(data_training$Alcohol, data_training$CAD, data_training$CHF)
interaction.plot(data_training$Alcohol, data_training$Diabetes, data_training$CHF)


par(mfrow=c(1,1))
interaction.plot(data_training$Alcohol, data_training$CAD, data_training$CHF)

#Full model with interactions
mod_glm = glm(CHF~(Gender+Age+Race+BP+Chol+Diabetes+CAD+PA+BMI+Alcohol)^2,data = data_training,family = binomial(link = "logit"))
summary(mod_glm)

mod_step = step(mod_glm,direction = "both")

mod_fin = glm(CHF ~ Age + Race + BP + Diabetes + CAD + PA + BMI + Alcohol + 
                Age:Diabetes + Age:CAD + Race:CAD + Race:Alcohol + Diabetes:CAD + 
                Diabetes:BMI + Diabetes:Alcohol + CAD:BMI + CAD:Alcohol + 
                PA:Alcohol + BMI:Alcohol,data=data_training,family = binomial(link = "logit"))
summary(mod_fin)



plot(mod_fin)


predict_glm = ifelse(predict(mod_fin,type = "response")>0.1,1,0)
conf_mat = table(data_training$CHF,predict_glm)
conf_mat


library(pROC)
library(ROCR)

predict_glm_test_1 = ifelse(predict(mod_fin,newdata = data_test, type = "response")>0.5,1,0)
predict_glm_test_2 = ifelse(predict(mod_fin,newdata = data_test, type = "response")>0.3,1,0)
predict_glm_test_3 = ifelse(predict(mod_fin,newdata = data_test, type = "response")>0.1,1,0)

conf_mat_test_1 = table(data_test$CHF,predict_glm_test_1)
conf_mat_test_2 = table(data_test$CHF,predict_glm_test_2)
conf_mat_test_3 = table(data_test$CHF,predict_glm_test_3)

conf_mat_test_1
conf_mat_test_2
conf_mat_test_3

plot(performance(prediction(predict_glm_test,data_test$CHF),"tpr","fpr"),main="ROCR")

auc_val_1 = auc(roc(data_test$CHF,predict_glm_test_1))
auc_val_1
auc_val_2 = auc(roc(data_test$CHF,predict_glm_test_2))
auc_val_2
auc_val_3 = auc(roc(data_test$CHF,predict_glm_test_3))
auc_val_3

par(mfrow=c(1,1))
plot(seq(0.1,0.5,0.2),c(auc_val_3,auc_val_2,auc_val_1),type = "line",col="red",xlab = "AUC value",ylab = "Cut-off value")

par(mfrow=c(1,2))
plot(seq(0.1,0.5,0.2),c(auc_val_3,auc_val_2,auc_val_1),type = "line",col="red",xlab = "AUC value",ylab = "Cut-off value")
plot(performance(prediction(predict_glm_test,data_test$CHF),"tpr","fpr"),main="ROCR")

library(dplyr)
odds_ratio = coef(mod_fin) %>% exp
odds_ratio


#Deviance test
dev_test_glm = mod_fin$null.deviance-mod_fin$deviance
p_val_dev_glm = 1 - pchisq(q=dev_test_glm,df=1)
p_val_dev_glm
#Reject null, thus, at least one covariate is related to the response


#Hosmer-Lemeshow goodness of fit test
install.packages("ResourceSelection")
library(ResourceSelection)
hoslem.test(data_training$CHF,fitted(mod_fin))
#Fail to reject the null, thus the model fits the data well


par(mfrow=c(1,2))
plot(predict(mod_fin),residuals(mod_fin),xlab = "Fitted values",ylab = "Residuals",main = "Residuals vs Fitted values")
abline(0,0,col="red")
qqnorm(residuals.glm(mod_fin))
qqline(residuals.glm(mod_fin),col="red")
