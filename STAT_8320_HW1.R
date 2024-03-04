library(dplyr)
library(matlib)
install.packages("pracma")
library(pracma)
install.packages("numDeriv")
library(numDeriv)
data_4 = read.table("C:\\Users\\Jayaditya Nath\\Downloads\\S24hw1pr4.txt")
data_4 = as.data.frame(data_4 %>% rename("X"="V1","Y"="V2"))
data_4
typeof(data_4$X)
model = function(theta_vec,x)
{
  theta_0 = theta_vec[1]
  theta_1 = theta_vec[2]
  theta_2 = theta_vec[3]
  y = (theta_0 + theta_1*x)/(1 + theta_2*exp(0.4*x))
  return(y)
}

initial_theta_curl = c(14,-6,4)
plot(data_4$X,data_4$Y)
lines(data_4$X,model(initial_theta_curl,data_4$X))


#Gauss_Newton Algorithm

Y = as.matrix(data_4$Y,nrow=250)
theta_curl = initial_theta_curl

for(iter in 1:10)
{
  f_theta = as.matrix(model(theta_curl,data_4$X),nrow=250)
  col_1 = eval(D(expression((theta_0 + theta_1*x)/(1 + theta_2*exp(0.4*x))),'theta_0'),list(theta_2=theta_curl[3],x=data_4$X))
  col_2 = eval(D(expression((theta_0 + theta_1*x)/(1 + theta_2*exp(0.4*x))),'theta_1'),list(theta_2=theta_curl[3],x=data_4$X))
  col_3 = eval(D(expression((theta_0 + theta_1*x)/(1 + theta_2*exp(0.4*x))),'theta_2'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x=data_4$X))
  F_theta = cbind(col_1,col_2,col_3)
  
  delta = pinv(t(F_theta) %*% F_theta) %*% t(F_theta) %*% (Y - f_theta)
  theta_curl = as.matrix(theta_curl) + delta
  
  if(max(abs(delta))<0.0001)
  {
    break
  }
}

print(iter)
print(paste0("The final estimates of theta are : ",theta_curl[1],",",theta_curl[2]," and ",theta_curl[3]," respectively."))
print(paste0("The algorithm converges at the ",iter," th iteration."))
print(paste0("The value of sigma_squared is : ",(t(Y - f_theta) %*% (Y - f_theta))/(nrow(data_4)-dim(theta_curl)[1])))

print("The variance-covariance matrix is of the form : ")
as.numeric((t(Y - f_theta) %*% (Y - f_theta))/(nrow(data_4)-dim(theta_curl)[1]))*as.matrix(solve(crossprod(F_theta)))
print(paste0("The value of the objective function at convergence is : ",sum((Y-f_theta)^2)))
Y-f_theta
par(mfrow=c(1,2)) 
plot(data_4$X,data_4$Y,main = "Curve fitted based on guess estimates")
lines(data_4$X,model(initial_theta_curl,data_4$X),col="green",lwd=2)
plot(data_4$X,data_4$Y,main="Curve fitted based on final estimates")
lines(data_4$X,model(theta_curl,data_4$X),col="magenta",lwd=2)












data_3_ini = read.table("C:\\Users\\Jayaditya Nath\\Downloads\\S24hw1pr3.txt",sep = ",")
data_3 = as.data.frame(data_3_ini[2:21,] %>% rename("X"="V1","Y"="V2"))
data_3 
typeof(data_3$X)
model_3 = function(theta_vec,x)
{
  theta_0 = theta_vec[1]
  theta_1 = theta_vec[2]
  theta_2 = theta_vec[3]
  y = (theta_0)/(1 + exp(-theta_1*(x-theta_2)))
  return(y)
}

theta_0_candidate = c(7.30, 7.41, 7.52, 7.63, 7.74, 7.86, 7.97, 8.08, 8.19, 8.30)
theta_1_candidate = c(0.30, 0.34, 0.39, 0.43, 0.48, 0.52, 0.57, 0.61, 0.66, 0.70)
theta_2_candidate = c(0.80, 0.84, 0.89, 0.93, 0.98, 1.02, 1.07, 1.11, 1.16, 1.20)

best_est_sse = 1000000
for (theta_0 in theta_0_candidate) {
  for (theta_1 in theta_1_candidate) {
    for (theta_2 in theta_2_candidate) {
      sse_est = sum((as.double(data_3$Y) - model_3(c(theta_0,theta_1,theta_2),as.double(data_3$X)))^2)
      
      if (sse_est < best_est_sse) {
        best_est_sse <- sse_est
        theta_curl = c(theta_0, theta_1, theta_2)
      }
    }
  }
}

print(paste0("The best estimates of theta are : ",theta_curl[1],",",theta_curl[2]," and ",theta_curl[3]," respectively."))
print(paste0("The SSE of the best estimates of theta is : ",best_est_sse))

print("The derivates of the mean function with respect to theta_0, theta_1 and theta_2 are : ")
D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0')
D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1')
D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2')

#Hessian check for part (e) 

hess_mat = matrix(rep(0,9),nrow = 3)

for(i in 1:dim(data_3)[1])
{
  dat <<- as.double(data_3$X)[i]
  new_mod = function(theta_new)
  {
    return(theta_new[1])/(1 + exp(-theta_new[2]*(dat-theta_new[3])))
    
  }
  new_mod(theta_curl)
  
  hess_mat = hess_mat + ((as.double(data_3$Y)[i] - model_3(theta_curl,as.double(data_3$X)[i])) * (hessian(func = new_mod,c(7.74 ,0.48, 0.98))))
}


gradient_new_rap = matrix(c(sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_0'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
),sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_1'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
),sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_2'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
)),nrow=3)
col_1_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
col_2_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
col_3_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
F_theta_new_rap = cbind(col_1_new_rap,col_2_new_rap,col_3_new_rap)
hess_mat_final = 2 * ((t(F_theta_new_rap)%*%F_theta_new_rap)-hess_mat)
alpha_learn_new_rap = 1
theta_curl_est_new_rap = as.matrix(theta_curl,nrow=3) - (alpha_learn_new_rap * inv(hess_mat_final) %*% gradient_new_rap)
print(paste0("The final estimates are : ",theta_curl_est_new_rap[1],", ",theta_curl_est_new_rap[2]," and ",theta_curl_est_new_rap[3]," respectively."))






Y = as.matrix(as.double(data_3$Y),nrow=20)
theta_curl_gauss_new = theta_curl

f_theta_gauss_new = as.matrix(model_3(theta_curl_gauss_new,as.double(data_3$X)),nrow=20)
col_1_gauss_new = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_curl_gauss_new[2],theta_2=theta_curl_gauss_new[3],x=as.double(data_3$X)))
col_2_gauss_new = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_curl_gauss_new[2],theta_2=theta_curl_gauss_new[3],x=as.double(data_3$X)))
col_3_gauss_new = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_curl_gauss_new[1],theta_1=theta_curl_gauss_new[2],theta_2=theta_curl_gauss_new[3],x=as.double(data_3$X)))
F_theta_gauss_new = cbind(col_1_gauss_new,col_2_gauss_new,col_3_gauss_new)
  
delta_gauss_new = pinv(t(F_theta_gauss_new) %*% F_theta_gauss_new) %*% t(F_theta_gauss_new) %*% (Y - f_theta_gauss_new)
theta_curl_gauss_new = as.matrix(theta_curl_gauss_new) + delta_gauss_new
  
if(max(abs(delta_gauss_new))<0.0001)
  {
    break
  }

print(paste0("The final estimates of theta are : ",theta_curl_gauss_new[1],",",theta_curl_gauss_new[2]," and ",theta_curl_gauss_new[3]," respectively."))
print(paste0("The Jacobian matrix is : ",F_theta_gauss_new))
print(paste0("The function value vector is : ",f_theta_gauss_new))
print(paste0("The inverse of cross-product of the Jacobian matrix is : ",pinv(t(F_theta_gauss_new) %*% F_theta_gauss_new)))
print(paste0("The increment vector is : ",delta_gauss_new))

theta_curl_gauss_new
f_theta_gauss_new
F_theta_gauss_new
pinv(t(F_theta_gauss_new) %*% F_theta_gauss_new)
delta_gauss_new


eval(D(expression(y-(theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(x=as.double(data_3$X),theta_1=theta_curl[2],theta_2=theta_curl[3]))
eval(D(expression(y-(theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3]))
eval(D(expression(y-(theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3]))

eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_0'),list(x=as.double(data_3$X),theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_1'),list(x=as.double(data_3$X),theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_2'),list(x=as.double(data_3$X),theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))

D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_1')

gradient = matrix(c(sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_0'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
),sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_1'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
),sum(eval(D(expression((y-(theta_0)/(1 + exp(-theta_1*(x-theta_2))))^2),'theta_2'),list(x=as.double(data_3$X),theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],y=as.double(data_3$Y)))
)),nrow=3)
alpha_learn = 1
theta_curl_est = as.matrix(theta_curl,nrow=3) - (alpha_learn*(diag(1,nrow=length(theta_curl),ncol=3) %*% gradient))
print(paste0("The final estimates are : ",theta_curl_est[1],", ",theta_curl_est[2]," and ",theta_curl_est[3]," respectively."))


col_1_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
col_2_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
col_3_new_rap = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_curl[1],theta_1=theta_curl[2],theta_2=theta_curl[3],x=as.double(data_3$X)))
F_theta_new_rap = cbind(col_1_new_rap,col_2_new_rap,col_3_new_rap)
2 * t(F_theta_new_rap) %*% F_theta_new_rap






theta_est_conv = c(7.6, 0.5, 1)

col_1_ci = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=as.double(data_3$X)))
col_2_ci = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=as.double(data_3$X)))
col_3_ci = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_est_conv[1],theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=as.double(data_3$X)))
F_theta_ci = cbind(col_1_ci,col_2_ci,col_3_ci)
f_theta_ci = as.matrix(model_3(theta_est_conv,as.double(data_3$X)),nrow=20)
sigma = sqrt((t(Y - f_theta_ci) %*% (Y - f_theta_ci))/(nrow(data_3)-length(theta_est_conv)))
sigma
print(paste0("95% confidence interval for theta_1 is (",theta_est_conv[2]-(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_ci) %*% F_theta_ci)[2,2])),",",theta_est_conv[2]+(qt(0.025,17,lower.tail = F)*sigma[1,1]*sqrt(inv(t(F_theta_ci) %*% F_theta_ci)[2,2])),")"))




f_theta_ci_2 = model_3(theta_est_conv,2.7)
col_1_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
col_2_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
col_3_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_est_conv[1],theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
F_theta_ci_2 = cbind(col_1_ci_2,col_2_ci_2,col_3_ci_2)
sigma_2 = sqrt((t(Y - f_theta_ci_2) %*% (Y - f_theta_ci_2))/(nrow(data_3)-length(theta_est_conv)))
print(paste0("95% confidence interval for Y|X=2.7 is (",f_theta_ci_2-(qt(0.025,17,lower.tail = F)*sigma_2[1,1]*sqrt(F_theta_ci_2 %*% pinv(t(F_theta_ci) %*% F_theta_ci) %*% t(F_theta_ci_2))),",",f_theta_ci_2+(qt(0.025,17,lower.tail = F)*sigma_2[1,1]*sqrt(F_theta_ci_2 %*% pinv(t(F_theta_ci_2) %*% F_theta_ci_2) %*% t(F_theta_ci_2))),")"))

print(paste0("95% prediction interval for Y|X=2.7 is (",f_theta_ci_2-(qt(0.025,17,lower.tail = F)*sigma_2[1,1]*sqrt(1 + (F_theta_ci_2 %*% pinv(t(F_theta_ci) %*% F_theta_ci) %*% t(F_theta_ci_2)))),",",f_theta_ci_2+(qt(0.025,17,lower.tail = F)*sigma_2[1,1]*sqrt(1 + (F_theta_ci_2 %*% pinv(t(F_theta_ci) %*% F_theta_ci) %*% t(F_theta_ci_2)))),")"))


h_theta_given = theta_est_conv[3]/theta_est_conv[2]
col_1_test = eval(D(expression(theta2/theta1),'theta0'),list(theta_est_conv[1],theta_est_conv[2],theta_est_conv[3]))
col_2_test = eval(D(expression(theta2/theta1),'theta1'),list(theta_est_conv[1],theta_est_conv[2],theta_est_conv[3]))
col_3_test = eval(D(expression(theta2/theta1),'theta2'),list(theta_est_conv[1],theta_est_conv[2],theta_est_conv[3]))
h_vect = c(col_1_test,col_2_test,col_3_test)
sigma = sqrt((t(Y - f_theta_ci) %*% (Y - f_theta_ci))/(nrow(data_3)-length(theta_est_conv)))
col_1_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_0'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
col_2_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_1'),list(theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
col_3_ci_2 = eval(D(expression((theta_0)/(1 + exp(-theta_1*(x-theta_2)))),'theta_2'),list(theta_0=theta_est_conv[1],theta_1=theta_est_conv[2],theta_2=theta_est_conv[3],x=2.7))
F_theta_ci_2 = cbind(col_1_ci_2,col_2_ci_2,col_3_ci_2)

print(paste0("The confidence interval for the given parameter combination is :(",h_theta_given - (qt(0.025,17,lower.tail = F)*sigma*sqrt(t(h_vect) %*% pinv(t(F_theta_ci_2) %*% F_theta_ci_2) %*% h_vect))
 ,",",h_theta_hat + (qt(0.025,17,lower.tail = F)*sigma*sqrt(t(h_vect) %*% pinv(t(F_theta_ci_2) %*% F_theta_ci_2) %*% h_vect))
,")"))



R_sq_pseudo = 1 - (sum((as.double(data_3$Y) - (as.matrix(model_3(theta_est_conv,as.double(data_3$X)),nrow=20)))^2)/sum((as.double(data_3$Y) - mean(as.double(data_3$Y)))^2))
print(paste0("The value of pseudo R_squared is : ",R_sq_pseudo))


F_theta_hat = cbind(col_1_ci,col_2_ci,col_3_ci)

mat_hat = F_theta_hat %*% inv(t(F_theta_hat) %*% F_theta_hat) %*% t(F_theta_hat)

for(i in seq_along(diag(mat_hat)))
{
  if(mat_hat[i,i]>((2*3)/20))
  {
    print(paste0("The ",i," th observation exhibits substantial leverage."))
  }
  else
  {
    print(paste0("The ",i," th observation does not exhibit substantial leverage."))
  }
}
cat("Hi",2)




