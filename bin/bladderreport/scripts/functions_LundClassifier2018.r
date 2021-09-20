library(ranger)

# load the classifiers
# load("./bin/scripts/classifiers_LundClassifier2018.Rdata")

# function to calculate the scores form SwitchBox classifiers
sample_score <- function(SampleData, classifier){
  trues <- SampleData[classifier$TSPs[,1],]  > SampleData[classifier$TSPs[,2],]
  score <- sum(classifier$score[trues])/sum(classifier$score)
  return(score)
}

# function to remove the rules from SwitchBox classifier when the genes that involved in the rules are not in the sample
keep_rules <- function(SampleData, classifier){
  # each rule has two genes, both of them should be in the sample data
  # otherwise this rule will be remove to avoid any errors
  sum  <- classifier$TSPs[,1] %in% rownames(SampleData) + classifier$TSPs[,2] %in% rownames(SampleData)
  
  # subset the classifier's rules
  classifier$name  <- paste0("TSPs", sum(sum==2))
  classifier$score <- classifier$score[sum==2]
  classifier$TSPs  <- classifier$TSPs[sum==2,]
  return(classifier)
}

# Random forest classifier prediction
# function to get the prediction and the scores from the Random forest classifier
RF_classify <- function(SampleData){
  
  require(ranger)
  
  # check the input data format # the input data will be converted to a data frame
  if(is.data.frame(SampleData)){    
    CurData <- SampleData[RF_InputGenes[,1],] > SampleData[RF_InputGenes[,2], , drop=FALSE]
    rownames(CurData) <- RF_InputGenes[,3] 
  }
  
  if(is.vector(SampleData)){ 
    CurData <- as.data.frame(SampleData[RF_InputGenes[,1]] > SampleData[RF_InputGenes[,2]], stringsAsFactors=FALSE)
    rownames(CurData) <- RF_InputGenes[,3]
  }
  
  
  if(is.matrix(SampleData)){  
    CurData <- as.data.frame(SampleData[RF_InputGenes[,1],] > SampleData[RF_InputGenes[,2],], stringsAsFactors=FALSE)
    rownames(CurData) <- RF_InputGenes[,3]
  }
  
  if (!is.data.frame(SampleData) & !is.vector(SampleData) & !is.matrix(SampleData)) {
    print(class(SampleData))
    print("Bad Format")
    return(NULL)
  }
  
  # prediction function
  Prediction <- predict(RF_classifier, as.data.frame(t(CurData), stringsAsFactors=FALSE))$predictions
  
  # Get the prediction based on the highest scores 
  CurPred <- colnames(Prediction)[apply(Prediction,1,which.max)]
  
  # change the names to Swedish
  CurPred <- gsub("Mes","Mesenkymal-lik",
                  gsub("Basal","Basal/Skivepitel-lik",
                       gsub("ScNE","Neuroendokrin-lik",
                            gsub("GU","Genomiskt Instabil",
                                 gsub("Uro","Urotel-lik ",
                                      CurPred)))))
  
  # make dataframe with two columns
  CurPred <- cbind(CurPred,CurPred)
  
  # change the UroA, UroB, and UroC to Uro in the first column (the column of the main subtypes)
  CurPred[grepl("Urotel-lik", CurPred[,1]), 1] <- "Urotel-lik"
  colnames(CurPred) <- c("Main_subtype","Subtype")
  
  # Prepare the scores
  # change the names to Swedish
  colnames(Prediction) <- gsub("Mes","Mesenkymal-lik",
                               gsub("Basal","Basal/Skivepitel-lik",
                                    gsub("ScNE","Neuroendokrin-lik",
                                         gsub("GU","Genomiskt Instabil",
                                              gsub("Uro","Urotel-lik ",
                                                   colnames(Prediction))))))
  
  # returning the subtypes and the scores
  return(list(Subtype=as.data.frame(CurPred, stringsAsFactors=FALSE, drop=F),
              Scores=Prediction))
}


# function to get the scores for main subtypes using SwitchBox classifiers
SwitchBox_classify <- function (SampleData){
  
  # check the input format
  # it will be converted to dataframe
  if(is.vector(SampleData)){
    SampleData <- as.data.frame(SampleData, stringsAsFactors=FALSE)
  }
  
  if(is.matrix(SampleData)){
    SampleData <- as.data.frame(SampleData, stringsAsFactors=FALSE)
  }
  
  if(!is.data.frame(SampleData) & !is.matrix(SampleData) & !is.vector(SampleData) ){    
    print("Bad format")
    return(NULL)
  }
  
  # check the gene rules for the main subtypes classifiers
  classifier_Uro_temp  <- keep_rules(SampleData, classifier_Uro)
  classifier_GU_temp   <- keep_rules(SampleData, classifier_GU)
  classifier_Basal_temp<- keep_rules(SampleData, classifier_Basal)
  classifier_Mes_temp  <- keep_rules(SampleData, classifier_Mes)
  classifier_ScNE_temp <- keep_rules(SampleData, classifier_ScNE)
  
  # check the gene rules for the Uro subclasses classifiers
  classifier_UroA_temp <- keep_rules(SampleData, classifier_UroA) 
  classifier_UroB_temp <- keep_rules(SampleData, classifier_UroB) 
  classifier_UroC_temp <- keep_rules(SampleData, classifier_UroC) 
  
  # check the gene rules the Grade classifier
  classifier_Grade_temp <- keep_rules(SampleData, classifier_Grade)
  
  # run the classifiers to get the scores
  # if the input is for one sample only
  if(ncol(SampleData)==1){
    
    # run the main subtypes classifiers
    Switchbox_main_subtype        <- rep(NA, 5)
    names(Switchbox_main_subtype) <- c("Uro","GU","Basal","Mes","ScNE")
    Switchbox_main_subtype["Uro"]  <- sample_score(SampleData, classifier_Uro_temp)
    Switchbox_main_subtype["GU"]   <- sample_score(SampleData, classifier_GU_temp)
    Switchbox_main_subtype["Basal"]<- sample_score(SampleData, classifier_Basal_temp)
    Switchbox_main_subtype["Mes"]  <- sample_score(SampleData, classifier_Mes_temp)
    Switchbox_main_subtype["ScNE"] <- sample_score(SampleData, classifier_ScNE_temp)
    
    # run the Uro subclasses classifiers
    Switchbox_Uro_subclass         <- rep(NA, 3)
    names(Switchbox_Uro_subclass)  <- c("UroA","UroB","UroC")
    Switchbox_Uro_subclass["UroA"] <- sample_score(SampleData, classifier_UroA_temp) 
    Switchbox_Uro_subclass["UroB"] <- sample_score(SampleData, classifier_UroB_temp) 
    Switchbox_Uro_subclass["UroC"] <- sample_score(SampleData, classifier_UroC_temp) 
    
    # run the Grade classifier
    Switchbox_Grade   <- c("Grade"=NA)
    Switchbox_Grade[] <- sample_score(SampleData, classifier_Grade_temp)
    
    # assign a sample subtype based on the highest score
    subtype     <- names(Switchbox_main_subtype)[which.max(Switchbox_main_subtype)]
    subtype_Uro <- names(Switchbox_main_subtype)[which.max(Switchbox_main_subtype)]
    
    # if the sample has Uro as a main subtype then get the Uro subclass prediction
    if(subtype=="Uro"){
      # assign a sample Uro subclass based on the highest score among the (UroA, UroB, UroC classifiers)
      subtype_Uro <- names(Switchbox_Uro_subclass)[which.max(Switchbox_Uro_subclass)]
    }
    
    # assign a sample molecular gerade 
    # if score > 0.5 then Grade 3, if not then it is not Grade 3
    Grade       <- if (Switchbox_Grade > 0.5) {"grad 3"} else {"ej grad 3"}
    
    # change the names to Swedish
    subtype <- gsub("Mes","Mesenkymal-lik",
                    gsub("Basal","Basal/Skivepitel-lik",
                         gsub("ScNE","Neuroendokrin-lik",
                              gsub("GU","Genomiskt Instabil",
                                   gsub("Uro","Urotel-lik",
                                        subtype)))))
    
    subtype_Uro <- gsub("Mes","Mesenkymal-lik",
                        gsub("Basal","Basal/Skivepitel-lik",
                             gsub("ScNE","Neuroendokrin-lik",
                                  gsub("GU","Genomiskt Instabil",
                                       gsub("Uro","Urotel-lik ",
                                            subtype_Uro)))))
    
    # prepare the scores for the output
    scores <- c(Switchbox_main_subtype=Switchbox_main_subtype,
                Switchbox_Uro_subclass=Switchbox_Uro_subclass,
                Switchbox_Grade=Switchbox_Grade)
    
    names(scores) <- c("Urotel-lik","Genomiskt Instabil","Basal/Skivepitel-lik","Mesenkymal-lik","Neuroendokrin-lik",
                       "Urotel-lik A","Urotel-lik B","Urotel-lik C")
    
    
    return(list(Subtype=as.data.frame(cbind(Main_subtype=subtype, 
                                            Subtype=subtype_Uro), drop=FALSE, stringsAsFactors=FALSE),
                Scores=as.matrix(t(scores)), 
                Grade=Grade))
    
    # if the input has >1 sample then repeat the code that for the one sample
  } else if(ncol(SampleData)>1){
    
    temp <- lapply(1:ncol(SampleData),function(x){
      Switchbox_main_subtype        <- rep(NA, 5)
      names(Switchbox_main_subtype) <- c("Uro","GU","Basal","Mes","ScNE")
      Switchbox_main_subtype["Uro"]  <- sample_score(SampleData[,x,drop=FALSE], classifier_Uro_temp)
      Switchbox_main_subtype["GU"]   <- sample_score(SampleData[,x,drop=FALSE], classifier_GU_temp)
      Switchbox_main_subtype["Basal"]<- sample_score(SampleData[,x,drop=FALSE], classifier_Basal_temp)
      Switchbox_main_subtype["Mes"]  <- sample_score(SampleData[,x,drop=FALSE], classifier_Mes_temp)
      Switchbox_main_subtype["ScNE"] <- sample_score(SampleData[,x,drop=FALSE], classifier_ScNE_temp)
      
      #
      Switchbox_Uro_subclass         <- rep(NA, 3)
      names(Switchbox_Uro_subclass)  <- c("UroA","UroB","UroC")
      Switchbox_Uro_subclass["UroA"] <- sample_score(SampleData[,x,drop=FALSE], classifier_UroA_temp) 
      Switchbox_Uro_subclass["UroB"] <- sample_score(SampleData[,x,drop=FALSE], classifier_UroB_temp) 
      Switchbox_Uro_subclass["UroC"] <- sample_score(SampleData[,x,drop=FALSE], classifier_UroC_temp) 
      
      #
      Switchbox_Grade   <- c("Grade"=NA)
      Switchbox_Grade[] <- sample_score(SampleData[,x,drop=FALSE], classifier_Grade_temp)
      
      subtype     <- names(Switchbox_main_subtype)[which.max(Switchbox_main_subtype)]
      subtype_Uro <- names(Switchbox_main_subtype)[which.max(Switchbox_main_subtype)]
      
      if(subtype=="Uro"){
        subtype_Uro <- names(Switchbox_Uro_subclass)[which.max(Switchbox_Uro_subclass)]
      }
      
      Grade       <- if (Switchbox_Grade > 0.5) {"grad 3"} else {"ej grad 3"}
      
      subtype<-gsub("Mes","Mesenkymal-lik",
                    gsub("Basal","Basal/Skivepitel-lik",
                         gsub("ScNE","Neuroendokrin-lik",
                              gsub("GU","Genomiskt Instabil",
                                   gsub("Uro","Urotel-lik",
                                        subtype)))))
      
      subtype_Uro<-gsub("Mes","Mesenkymal-lik",
                        gsub("Basal","Basal/Skivepitel-lik",
                             gsub("ScNE","Neuroendokrin-lik",
                                  gsub("GU","Genomiskt Instabil",
                                       gsub("Uro","Urotel-lik ",
                                            subtype_Uro)))))
      
      return(list(Subtype=as.data.frame(cbind(Main_subtype=subtype, 
                                              Subtype=subtype_Uro),drop=FALSE,stringsAsFactors=FALSE),
                  Scores=c(Switchbox_main_subtype=Switchbox_main_subtype,
                           Switchbox_Uro_subclass=Switchbox_Uro_subclass,
                           Switchbox_Grade=Switchbox_Grade),
                  Grade=Grade))
    }
    )
    
    subtype <- do.call("rbind",lapply(temp,function(x){x[[1]]}))
    scores  <- do.call("rbind",lapply(temp,function(x){x[[2]]}))
    Grade   <- do.call("rbind",lapply(temp,function(x){x[[3]]}))
    
    colnames(scores)<-c("Urotel-lik",
                        "Genomiskt Instabil",
                        "Basal/Skivepitel-lik",
                        "Mesenkymal-lik",
                        "Neuroendokrin-lik",
                        "Urotel-lik A",
                        "Urotel-lik B",
                        "Urotel-lik C",
                        "Grade")
    
    return(list(Subtype=as.data.frame(subtype, drop=FALSE, stringsAsFactors=FALSE),
                Scores=scores,
                Grade=Grade,
                Grader_num=NA))
  } 
}


# Final combined classification
# This function use the SwitchBox predictions and the Random Forest predictions  
LundClassifier <- function(SampleData){
  
  # check the input format
  # it will convert it to dataframe
  if(is.vector(SampleData)){ 
    SampleData <- as.data.frame(SampleData, stringsAsFactors=FALSE)
  }
  
  if(is.matrix(SampleData)){  
    SampleData <- as.data.frame(SampleData, stringsAsFactors=FALSE)
  }  
  
  if(!is.data.frame(SampleData) & !is.matrix(SampleData) & !is.vector(SampleData) ){    
    print("Bad format")
    return(NULL)
  }
  
  # classify using Random forest classifier
  RF_pred <- RF_classify(SampleData)
  
  # classifiy using SwitchBox classifiers
  SB_pred <- SwitchBox_classify(SampleData)
  
  # report the Consensus of the SwitchBox and Random Forest models 
  # the prediction for the main subtype is reported under Main_subtype (the 5 main subtype)
  # the prediction for the subclass is reported under Subtype (4 main subtypes with Uro subclasses)
  # if the two models did not agree then "Ej Bedömbar" will be reported
  Consensus <- do.call("rbind",lapply(1:ncol(SampleData),function(x){
    
    LenL<-(length(unique(  c(RF_pred[["Subtype"]][x,"Main_subtype"],SB_pred[["Subtype"]][x,"Main_subtype"])  ) )!=1)
    LenS<-(length(unique(  c(RF_pred[["Subtype"]][x,"Subtype"],SB_pred[["Subtype"]][x,"Subtype"])  ) )!=1)
    
    return(c(Main_subtype=if(LenL){
      paste(unique(c("Ej Bedömbar:",RF_pred[["Subtype"]][x,"Main_subtype"],"|",SB_pred[["Subtype"]][x,"Main_subtype"])),
            collapse="")}
      else{
        paste(unique(c(RF_pred[["Subtype"]][x,"Main_subtype"],SB_pred[["Subtype"]][x,"Main_subtype"])),
              collapse="|")},
      Subtype=if(LenS){
        paste(unique(c("Ej Bedömbar:",RF_pred[["Subtype"]][x,"Subtype"],"|",SB_pred[["Subtype"]][x,"Subtype"])),
              collapse="")}
      else{
        paste(unique(c(RF_pred[["Subtype"]][x,"Subtype"],SB_pred[["Subtype"]][x,"Subtype"])),collapse="|")}, 
      Grade=as.character(SB_pred[["Grade"]][x])))
  }) )
  
  Report<-do.call("rbind",lapply(1:ncol(SampleData),function(x){
    
    LenL<-(length(unique(  c(RF_pred[["Subtype"]][x,"Main_subtype"],SB_pred[["Subtype"]][x,"Main_subtype"])  ) )!=1)  
    LenS<-(length(unique(  c(RF_pred[["Subtype"]][x,"Subtype"],SB_pred[["Subtype"]][x,"Subtype"])  ) )!=1) 
    
    return(c(Main_subtype=if(LenL){paste(unique(c("Ej Bedömbar:")),collapse="")}
             else{
               paste(unique(c(RF_pred[["Subtype"]][x,"Main_subtype"],SB_pred[["Subtype"]][x,"Main_subtype"])),collapse="|")},
             Subtype=if(LenS){paste(unique(c("Ej Bedömbar:")),collapse="")}
             else{
               paste(unique(c(RF_pred[["Subtype"]][x,"Subtype"],SB_pred[["Subtype"]][x,"Subtype"])),collapse="|")}, 
             Grade=as.character(SB_pred[["Grade"]][x])))
  }) )
  
  return(list(Report      =as.data.frame(Report, drop=FALSE, stringsAsFactors=FALSE),
              Consensus   =as.data.frame(Consensus, drop=FALSE, stringsAsFactors=FALSE),
              RandomForest=RF_pred,
              SwitchBox   =SB_pred))
}