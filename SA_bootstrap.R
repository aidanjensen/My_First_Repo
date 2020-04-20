#This file contains code to load the data from a .csv downloaded from the "SA Hub for BI Team" Google Sheet
#This file also contains SA_boostrap() which will run a bootstrap analysis on the SA_data object
#Finally, this file contains graph.wins() and graph.mean.rsq() which take the output of SA_bootstrap and plot them quickly
#This code was written by Aidan Jensen :)

require(readr)
require(ggplot2)
require(dplyr)
require(bestNormalize)

#Load the data into SA_data for analysis
SA_data <- read_csv("~/Downloads/SA Hub for BI Team - Data for Analysis in R.csv")

#View the data to verify
View(SA_data)

#Make BoxCox transformed columns

#do a boxcox analysis on each column
SA.boxcox <-function(df){
bc.pos <- boxcox(df$positive)
bc.neg <- boxcox(df$negative)
bc.neu <- boxcox(df$neutral)
#use the lambdas from each boxcox analysis to add a transformed column
return(df %>% mutate(bc.positive = positive ^ bc.pos$lambda) %>% mutate(bc.negative = negative ^ bc.neg$lambda) %>% 
  mutate(bc.neutral = neutral ^ bc.neu$lambda))
}

SA.logit <- function(df){
  df$logit.pos <- log((df$positive/(1-df$positive)))
  df$logit.neg <- log((df$negative/(1-df$negative)))
  df$logit.neu <- log((df$neutral/(1-df$neutral)))
  return(df)
}




lm#The bootstrap function
SA_bootstrap <- function(df = SA_data,iterations = 100, training.data.fraction = .75, metaseed = 42, loud = FALSE){
  require(bestNormalize)
  df2 <- SA.boxcox(df)
  df3 <- SA.logit(df2)
  require(dplyr)
  
  #initialize some variables for our results table
  
  #wins
  lm.wins <- 0
  bc.lm.wins <- 0
  logit.lm.wins <-0
  VADER.wins <- 0
  NC.wins <- 0
  Bucket.wins <- 0
  
  #sum of adjusted R squared values
  lm.sum.adj.Rsq <- 0
  bc.lm.sum.adj.Rsq <- 0
  logit.lm.sum.adj.Rsq <- 0
  VADER.sum.adj.Rsq <- 0
  NC.sum.adj.Rsq <- 0
  Bucket.sum.adj.Rsq <- 0
  
  #Put all of these initialized variables into a tibble
  results.db <- as_tibble(data.frame(lm.wins, bc.lm.wins, logit.lm.wins, VADER.wins, NC.wins, Bucket.wins, lm.sum.adj.Rsq, bc.lm.sum.adj.Rsq, logit.lm.sum.adj.Rsq, VADER.sum.adj.Rsq, NC.sum.adj.Rsq, Bucket.sum.adj.Rsq))
  
  #use a metaseed to make sure we get the same bunch of seeds for replicability (default metaseed = 42)
  set.seed(metaseed)
  #grab a handful of seeds so we can add maximum randomization to our bootstrap
  seeds <- sample(x = 1:100000000, size = iterations)
  
  #here is the actual bootstrap
  
  #For each seed
  for(i in 1:length(seeds)){
    
    #print the iteration and set the RNG seed
    if(loud){
      (print(i))}
    
    set.seed(seeds[i])
    
    #grab the indices of data to be used in training set
    ind <- sample(1:nrow(df2), size = round(nrow(df3) * training.data.fraction))
    
    #Define our training and testing data
    training_data <- df3[ind,]
    test_data <- df3[-ind,]
    
    #Use the training data to fit a linear model on HumanScore using predictors neutral, negative, and positive from the training set
    fit <- lm(HumanScore ~ neutral + negative + positive, data = training_data)
    fit.bc <- lm(HumanScore ~ bc.positive + bc.negative + bc.neutral, data = training_data)
    fit.logit <- lm(HumanScore ~ logit.pos + logit.neg + logit.neu, data = training_data)
    
    #Use the coefficients calculated to predict the HumanScore of the test_data
    test_data$LinearModel <- predict(fit, newdata = test_data)
    test_data$BoxCoxLM <- predict(fit.bc, newdata = test_data)
    test_data$LogitLM <- predict(fit.logit, newdata = test_data)
    
    #Generate a summary of the performance of each SA method
       #This summary will outline the ability of each SA method to predict HumanScore using linear modeling
       #Linear modeling should be OK because a small change in the method's score should coincide with a small change in predicted HumanScore
          #We are mostly interested in this summary for the adjusted R-squared value that is produced
    lmfit <- summary(lm(HumanScore~ LinearModel, data = test_data))
    bc.lmfit <- summary(lm(HumanScore ~ BoxCoxLM, data = test_data))
    logit.lmfit <- summary(lm(HumanScore ~ LogitLM, data = test_data))
    VADERfit <- summary(lm(HumanScore ~ VADER, data = test_data))
    NCfit <- summary(lm(HumanScore ~ AWS_NaiveContinuous, data = test_data))
    Bucketfit <- summary(lm(HumanScore ~ AWS_Bucket, data = test_data))
    
    #create an object, v, which stores the adjusted R-squared values from each summary
    v <- c(lmfit$adj.r.squared, bc.lmfit$adj.r.squared, logit.lmfit$adj.r.squared, VADERfit$adj.r.squared, NCfit$adj.r.squared, Bucketfit$adj.r.squared)
    
    #Which method has the highest adjusted R squared for this iteration
       #This method is chosen as the "winner" of this iteration
    max.Rsq <- which.max(v)
    
    #Increment the winner's win count
    results.db[,max.Rsq] <- results.db[,max.Rsq] + 1
    
    #Add the sum of adjusted R-squared for each method
    results.db$lm.sum.adj.Rsq <- results.db$lm.sum.adj.Rsq + lmfit$adj.r.squared
    results.db$VADER.sum.adj.Rsq <- results.db$VADER.sum.adj.Rsq + VADERfit$adj.r.squared
    results.db$NC.sum.adj.Rsq <- results.db$NC.sum.adj.Rsq + NCfit$adj.r.squared
    results.db$Bucket.sum.adj.Rsq <- results.db$Bucket.sum.adj.Rsq + Bucketfit$adj.r.squared
    results.db$bc.lm.sum.adj.Rsq <- results.db$bc.lm.sum.adj.Rsq + bc.lmfit$adj.r.squared
    results.db$logit.lm.sum.adj.Rsq <- results.db$logit.lm.sum.adj.Rsq + logit.lmfit$adj.r.squared
    
    #end of iteration
  }
  #When all iterations of bootstrap completed, return the finished results.db containing wins and sum of adjusted R-squared
  return(results.db)
  
}



#The graph.wins function needs a results.db as an argument!
graph.wins <- function(boot.res.db){
  require(ggplot2)
  #let's tweak the column names a little bit
colnames(boot.res.db) <- c("LinearModel_Wins", "BoxCoxLM_Wins", "VADER_Wins", "NaiveContinuous_Wins", "Bucket_Wins",
                                              "LinearModel_MEAN.Adj.R^2","BoxCoxLM_MEAN.Adj.R^2", "VADER_MEAN.Adj.R^2", "NaiveContinuous_MEAN.Adj.R^2", "Bucket_MEAN.Adj.R^2")
 #label the methods
methods <- c("LinearModel","BoxCoxLM", "VADER", "NaiveContinuous", "Bucket")
 #transfer number of wins
wins <- c(boot.res.db$LinearModel_Wins,boot.res.db$BoxCoxLM_Wins, boot.res.db$VADER_Wins, boot.res.db$NaiveContinuous_Wins, boot.res.db$Bucket_Wins)

#format the data
x <- as_tibble(methods)
x <- cbind(x, wins)

#produce the barplot
ggplot(x, aes(x = reorder(methods, -wins) , y = wins))+geom_bar(stat = "identity")+ylab("# of iterations containing best adj. R^2")+xlab("Method")+ggtitle("Number of iteration wins per Sentiment Analysis method")+coord_flip()
}


graph.sum.rsq <- function(boot.res.db){
  require(ggplot2)
  #same process as previous function
  methods <- c("LinearModel","BoxCoxLM", "VADER", "NaiveContinuous", "Bucket")
  colnames(boot.res.db) <- c("LinearModel_Wins", "BoxCoxLM_Wins", "VADER_Wins", "NaiveContinuous_Wins", "Bucket_Wins",
                             "LinearModel_MEAN.Adj.R^2","BoxCoxLM_MEAN.Adj.R^2", "VADER_MEAN.Adj.R^2", "NaiveContinuous_MEAN.Adj.R^2", "Bucket_MEAN.Adj.R^2")  
  mean_rsq <- c(boot.res.db$`LinearModel_MEAN.Adj.R^2`, boot.res.db$`BoxCoxLM_MEAN.Adj.R^2`, boot.res.db$`VADER_MEAN.Adj.R^2`, boot.res.db$`NaiveContinuous_MEAN.Adj.R^2`, boot.res.db$`Bucket_MEAN.Adj.R^2`)
  x <- as_tibble(methods)
  x <- cbind(x, mean_rsq)
  ggplot(x, aes(x = reorder(methods, -mean_rsq) , y = mean_rsq))+geom_bar(stat = "identity")+ylab("Mean Adjusted R-Squared")+xlab("Method")+ggtitle("Mean Adjusted R-squared by Sentiment Analysis method")+coord_flip()

  }
