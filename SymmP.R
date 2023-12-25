library(lme4)
l11<-c()
l21<-c()
l31<-c()
in11<-c()
in21<-c()
in31<-c()
#--------------Compute prediction coverage probability and length of prediction set-----------

for (rep2 in c(1:100)){ #100 independent trials
  l1<-c()
  l2<-c()
  l3<-c()
  in1<-c()
  in2<-c()
  in3<-c()
for (rep in c(1:500)){ #500 test samples
  len1<-c()
  len2<-c()
  len3<-c()
  Indt1<-c()
  Indt2<-c()
  Indt3<-c()
  print(rep)
  data2<-matrix(sleepstudy$Reaction,nrow=10,ncol=18)
  samp<-c(1,sample(c(1:9),9,replace=FALSE)+1) #a permutation of days. We pick the first 6 entries to be the training data
  data1<-data2[samp,]
  a1<-sleepstudy$Reaction #sleepreaction time
  b1<-matrix(sleepstudy$Days,nrow=10,ncol=18) #number of days in the experiment
  b<-b1[samp,]
  
  
  X2<-c()
  #train data:
  for ( i in c(1:18)){
    X2<-c(X2,rep(data2[1,i],6)) #baseline time
  }
  X1<-c()
  for (i in c(1:18)){ #indices of traning data
    X1<-c(X1,b[2:7,i])
  }
  
  Y_train<-c() #response variable
  for (i in c(1:18)){
    Y_train<-c(Y_train,data1[2:7,i])
  }
  
  
  train<-data.frame(Y_train,X1,X2) #training the model
  fit<-lm(Y_train~.-1,data=train)  
  
  
  X2<-matrix(X2,nrow=6) #collection of training data
  X1<-matrix(X1,nrow=6)
  Y_train<-matrix(Y_train,nrow=6) 
  sd2<-c()
  coef<-matrix(nrow=18,ncol=2)
  res<-matrix(nrow=6,ncol=18)
  #----------Training the data--------------------------
  for (j in c(1:18)){
    R<-Y_train[,j]
    C1<-X1[,j]
    C2<-X2[,j]
    train<-data.frame(R,C1,C2)
    fit2<-lm(R~.-1,data=train)  
    sd2<-c(sd2,sd(R-fit2$coefficients[1]*C1-fit2$coefficients[2]*C2))
    res[,j]<-R-fit2$coefficients[1]*C1-fit2$coefficients[2]*C2
    su<-summary(fit2)
    coef[j,]<-fit2$coefficients
  }
  
  for ( i in which(coef[,1]<0)){
    coef[i,1]<-fit$coefficients[1]}
  #-------------X_1t,X_2t are calibration samples---
  X1_t<-matrix(nrow=3,ncol=18)
  for ( i in c(1:18)){
    X1_t[,i]<-b[8:10,i]
  }
  
  
  X2_t<-matrix(nrow=3,ncol=18)
  for ( i in c(1:18)){
    X2_t[,i]<-rep(data2[1,i],3)
  }
  
  
  
  Y_test<-matrix(nrow=3,ncol=18)
  for ( i in c(1:18)){
    Y_test[,i]<-data1[8:10,i]
  }
  
  
  #--------Predict---------
  Y_pred<-matrix(nrow=3,ncol=18)
  
  for (i in c(1:18)){
    Y_pred[,i]<-coef[i,1]*X1_t[,i]+coef[i,2]*X2_t[,i] #\hat \mu_k, j=1\cdots, 18 in our method
  }
  Y_pred2<-matrix(nrow=3,ncol=18)
  for (i in c(1:18)){
    Y_pred2[,i]<-fit$coefficients[1]*X1_t[,i]+fit$coefficients[2]*X2_t[,i]  #\hat \mu in conformal prediction
  }
  #Y_pred2
  seq0=seq(-200,500,length=5000) #test sample
  ind1<-c()
  ind2<-c()
  Y_test1<-Y_test
  c=2
  band<-matrix(nrow=3,ncol=18)
  
  for(i in c(1:3)){  #calculate the confidence band for linear regression
    for (j in c(1:18)){
      cov<-cbind(X1[,i],X2[,i])
      s<-t(cov)%*%cov 
      #solve(s)
      tz<-cbind(X1_t[i,j],X2_t[i,j])
      ma<-tz%*%solve(s)%*%t(tz)
      band[i,j]<-c*sqrt(ma)*sd2[j]
    }
  }
  #------------Random select one test point among 54 calibration sets and build prediction interval---------------------
  i=sample(c(1:3),size=1)
  j=sample(c(1:18),size=1)
      Y_test1<-Y_test
      true=Y_test[i,j]
      ind1<-c()
      ind2<-c()
      ind3<-c()
      for (r in c(1:length(seq0))){
        Y_test1[i,j]<-seq0[r]
        Mat<-Y_test1-Y_pred2*I(abs(Y_pred-Y_pred2)<band)-Y_pred*I(abs(Y_pred-Y_pred2)>band)
        Mat0<-Y_test1-Y_pred2  
  #----Calibration Estimators----------------
        Mat1<-Mat
        for (k in c(1:18)){
          Mat1[,k]<-abs(Mat[,k]-mean(Mat[,k]))/max(sd2[k],20)
          Mat0[,k]<-abs(Mat0[,k])
        }
        Mat3<-Mat0[1,setdiff(c(1:18),j)]
        Mat3<-c(Mat3,Mat0[i,j])
        if (sum(Mat1[i,j]<Mat1)<5){
          ind1<-c(ind1,1)
        }else{ind1<-c(ind1,0)}
        if (sum(Mat0[i,j]<Mat0)<5){
          ind2<-c(ind2,1)
        }else{ind2<-c(ind2,0)}
        if (sum(Mat3[length(Mat3)]<Mat3)<1){
          ind3<-c(ind3,1)
        }else{ind3<-c(ind3,0)}
      }
      #--------------------Record the length of the prediction set and the indicator of coverage------------
      len1<-c(len1,seq0[which(ind1==0)[length(which(ind1==0))]]-seq0[which(ind1==0)[1]]) #our length
      len2<-c(len2,seq0[which(ind2==0)[length(which(ind2==0))]]-seq0[which(ind2==0)[1]]) #conformal prediction
      len3<-c(len3,seq0[which(ind3==0)[length(which(ind3==0))]]-seq0[which(ind3==0)[1]]) #sub-sampling once
      if ((seq0[which(ind1==0)[1]]<=true)&(true<=seq0[which(ind1==0)[length(which(ind1==0))]])){ #our indicator on whether the prediction interval covers the true value
        Indt1<-c(Indt1,1)
      }else{
        Indt1<-c(Indt1,0)
      }
      if ((seq0[which(ind2==0)[1]]<=true)&(true<=seq0[which(ind2==0)[length(which(ind2==0))]])){ #indicator of conformal prediction
        Indt2<-c(Indt2,1)
      }else{Indt2<-c(Indt2,0)
      }
      if ((seq0[which(ind3==0)[1]]<=true)&(true<=seq0[which(ind3==0)[length(which(ind3==0))]])){ #indicator of sub-sampling
        Indt3<-c(Indt3,1)
      }else{
        Indt3<-c(Indt3,0)
      }
    #-----------Record the coverage probabilities and lengths of the prediction sets of 500 test samples--------------
      in1<-c(in1,Indt1)
      write.csv(in1,"in1.csv")  #Coverage: our method
      in2<-c(in2,Indt2)
      write.csv(in2,"in2.csv") #Coverage: conformal prediction
      in3<-c(in3,Indt3)
      write.csv(in3,"in3.csv") #Coverage: subsampling
      l1<-c(l1,len1)
      write.csv(l1,"l1.csv") #Length: our method
      l2<-c(l2,len2)
      write.csv(l2,"l2.csv") #Length: conformal prediction
      l3<-c(l3,len3)
      write.csv(l3,"l3.csv") #Length: sub-sampling
}
  #----------_Record the averaged lengths and stds of 100 trials------------------
  in11<-c(in11,mean(in1))  #averaged coverage: our method 
  write.csv(in11,"in11.csv")
  in21<-c(in21,mean(in2)) #averaged coverage: conformal prediction
  write.csv(in21,"in21.csv")
  in31<-c(in31,mean(in3)) #averaged coverage: sub-sampling
  write.csv(in31,"in31.csv")
  l11<-c(l11,mean(l1))   #averaged length: our method 
  write.csv(l11,"l11.csv")
  l21<-c(l21,mean(l2))   #averaged length: conformal prediction
  write.csv(l21,"l21.csv")
  l31<-c(l31,mean(l3))   #averaged length: sub-sampling
  write.csv(l31,"l31.csv")
}