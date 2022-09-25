
# HEART ATTACKS ANALYSIS AND PREDICTION

# 1 LIBRARIES

library( GGally)
library(MASS)
library( rms )
library(arm)
library(ResourceSelection)
library(pROC)
library(PRROC)

# 2 DATASET

dati = read.csv( "heart.csv", head = TRUE )

# We remove na values.

na.omit(dati)

# The dataset presents the following covariates:

# Age : Age of the patient
# Sex : Sex of the patient 1 = male; 0 = female
# exang: exercise induced angina (1 = yes; 0 = no)
# caa: number of major vessels (0-3)
# cp : Chest Pain type chest pain type
# Value 1: typical angina
# Value 2: atypical angina
# Value 3: non-anginal pain
# Value 4: asymptomatic
# trtbps : resting blood pressure (in mm Hg)
# chol : cholestoral in mg/dl fetched via BMI sensor
# fbs : (fasting blood sugar > 120 mg/dl) (1 = true; 0 = false)
# rest_ecg : resting electrocardiographic results
# Value 0: normal
# Value 1: having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV)
# Value 2: showing probable or definite left ventricular hypertrophy by Estes' criteria
# thalach : maximum heart rate achieved
# output : 0 (no heart disease) or 1 (heart disease)

# Given the fact that some features are categorical, we use the factor operator in order to work with them later.

fattoriali = c('sex', 'exng', 'cp', 'fbs', 'restecg')
dati[fattoriali] = lapply(dati[fattoriali], factor)

# We can use the package 'GGally'. We can easily vizualize the relationship of couple of variables, 
# their sample correlation and their approximated density function.

ggpairs(data = dati, title ="Relationships between predictors & response",
        lower = list(continuous=wrap("points", alpha = 0.5, size=0.1)))

# 3 FEATURE SELECTION

# We immediately remove the following variables because they seem difficult to understand: we tend to
# prefer a simple model (with a limited number of features).

dati = subset(dati, select = -c(caa,thall,oldpeak,slp))

# Most of the information provided by 'exng' is also provided by 'cp' (Chest pain). So, we remove 'exng' to prefer 
# model parsimony and to avoid data redundancy.

dati = subset(dati, select = -exng)

# We decide to visualize the output distribution dividing our dataset according to the sex of the patient.

par(mfrow = c(1,2))
table_output_women=table(dati$output[which(dati$sex==0)])
table_output_men=table(dati$output[which(dati$sex==1)])

barplot(table_output_women,names.arg=c('0-no stroke','1-stroke'), col = 'plum', xlab = 'Esito', ylab = 'Numerosità', main = "Distribuzione dell'output nelle donne")
barplot(table_output_men,names.arg=c('0-no stroke','1-stroke'), col = 'slategray1', xlab = 'Esito', ylab = 'Numerosità', main = "Distribuzione dell'output negli uomini")

stroke_women=table_output_women[2]/(table_output_women[1]+table_output_women[2])
stroke_men=table_output_men[2]/(table_output_men[1]+table_output_men[2])

# From our data, we notice women seem to be more exposed to heart diseases than men. Considering the structure of 
# our dataset, this may be due to the fact that most of the observations regard male patients. Therefore, even a small 
# number of heart diseases in the female category can have a huge impact on our model and our results.

# 4 REGRESSION MODEL

# We can fit the regression model and look at a summary of the estimated coefficients.

modello = glm( output ~ ., family = binomial( link = logit ), data=dati )
summary( modello )

# High p-values suggest that the covariates they refer to are not statistically significant. Therefore, we compute a stepwise 
# forward selection of variables removing the feature with the highest p.value (fbs ~ 0.81).

modello1 = update(modello, . ~ . -fbs)
summary( modello1 )

# We now remove 'restecg' because it seems not to be significant according to the Z-test computed by R.

modello2 = update(modello1, . ~ . -restecg)
summary( modello2 )

# Using the ANOVA test, we can contrast MODELLO and MODELLO2 in terms of fit and parsimony. THe variation of deviance from
# the models does not seem to be statistically significant (high p-value ~ 0.5332), therefore we prefer the second model
# which is simpler and more interpretable than MODELLO.

anova( modello2, modello, test = "Chisq" )

# We decide to mantain the covariate 'age' because, according to medical literature, it has a pretty huge role 
# in the onsets of heart diseases.

# We once again visualize the relationship of couple of variables and their correlation.

ggpairs(data = dati[,c('output', 'age', 'sex', 'cp', 'trtbps', 'chol', 'thalachh')], title ="Relationships between predictors & response",
        lower = list(continuous=wrap("points", alpha = 0.5, size=0.1)))

modello3 = update(modello2, . ~ . +age*trtbps)
summary( modello3 )

modello4 = update(modello2, . ~ . +age*chol)
summary( modello4 )

# We decide to introduce the variable 'age*thalachh' in order to weight the effect that age and thalachh have on the output.
# From the summary of MODELLO5 we have strong statistical evidence of the significance of this new variable.

modello5 = update(modello2, . ~ . +age*thalachh)
summary( modello5 )

# We also notice that the deviance and AIC value of MODELLO5 are smaller than the values of the previous model (MODELLO2). 

# Moreover, the high p-value (~ 0.758) obtained from the hoslem test confirms the fact that the model we are considering is 
# well performing. As a matter of fact, a high p-value suggests that the outputs computed by the model are not that different from the 
# empirical observed outputs. 

hoslem.test( modello5$y, fitted( modello5 ), g = 12) # p-value molto alto (0.758)

modello6 = update(modello5, . ~ . -chol)
summary( modello6 )

anova( modello6, modello5, test = "Chisq" )

modello7 = update(modello6, . ~ . -trtbps)

summary( modello7 )

anova( modello7, modello5, test = "Chisq" )

# MODELLO7 is the final model we obtained and we will analyze its coefficients and performance from now on.

# 5 COEFFICIENTS

# We compute the confidence regions for the coefficients of each covariate. As expected, with a confidence level of 95%,
# none of these intervals contains the value 0.

sm7=summary(modello7)
alpha=0.05
r_tot=length(coef(modello7))
for(i in (1:8)){
  ICmin=coef(modello7)[i]-coef(sm7)[i,2]*qnorm(1-alpha/2)
  ICmax=coef(modello7)[i]+coef(sm7)[i,2]*qnorm(1-alpha/2)
  IC=c(ICmin,ICmax)
  print(IC)}

# 6 PERFORMANCE

# We compute the sensitivity and specificity of our model in order to judge its performance and goodness of fit.

soglia = 0.5
valori.reali  = dati$output
valori.predetti = as.numeric( modello7$fitted.values > soglia ) # we convert boolean values into numeric ones
tab = table( valori.reali, valori.predetti )

# Sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# Specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita

# The obtained values seem pretty convincing even with a threshold equal to 0.5. We do not need to change it.

# 7 ROC CURVE

curva_roc <- roc.curve(scores.class0 = modello7$fitted, weights.class0=as.numeric(paste(dati$output)), curve=TRUE)
plot(curva_roc)

# 8 CROSS VALIDATION

# We implement a cross validation in order to validate our model. We divide our dataset into two different blocks:

# 1) the training set consists of 80% of the observations we have;

# 2) the test set consists of 20% of the observations we have.

# We repeat the procedure 20 times using the 'Brier Score' as metric to evaluate the performance of the model.

dim_tot = dim(dati)[1]

dim_train = ceiling(0.8*dim_tot)
dim_test = dim_tot-dim_train
err=rep(0,20)

for( i in (1:20)){
  train_rows=sample(seq(1,dim_tot),replace=FALSE,size=dim_train)
  train_df=dati[train_rows,]
  test_df=dati[-train_rows,]
  modello_temp=glm(output~ age + sex + cp+ thalachh+ thalachh*age, family=binomial(link=logit),data=train_df)
  valori_predetti=predict(modello_temp,test_df,se=TRUE)
  err[i]=mean((test_df$output - modello_temp$family$linkinv(valori_predetti$fit))^2)
}
max(err)
min(err)
mean(err)

# We get the following results:

# - max err = 0.187

# - min err = 0.098

# - average err = 0.145

# Therefore, we can state the model is well performing other than easily interpretable.












