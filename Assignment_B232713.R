############ Question (1) ############

# Install necessary packages and load them accordingly
install.packages("tidyverse")
library("tidyverse")
library("readxl")
library("dplyr")
library("ggplot2")
library("psych")

# Download the excel files from Github and read it in R Studio
download.file("https://github.com/s2465202/Probability-Statistics_Assignment/raw/main/biomarkers.xlsx", destfile = "biomarkers.xlsx",mode = "wb")
biomarker <- read_excel("biomarkers.xlsx")
download.file("https://github.com/s2465202/Probability-Statistics_Assignment/raw/main/covariates.xlsx", destfile = "covariates.xlsx",mode = "wb")
covariates <- read_excel("covariates.xlsx")

# Create new columns to separate the patient ID and time of blood sampling and convert the table from long to wide format in order to have 1 row per patient
biomarker_name <- colnames(biomarker[,2:10])
biomarker <- biomarker %>%
  separate_wider_delim(Biomarker,delim="-",names=c("PatientID","SampleTime")) %>%
  pivot_wider(names_from=SampleTime,values_from=all_of(biomarker_name),names_vary="slowest")

# Review data structure and data types for biomarker and covariates to see if any conversion if required before joining up the two datasets and using them for analysis
str(biomarker)
str(covariates)

# All of the columns in the two datasets are numerical except for biomarker$PatientID. Convert biomarker$PatientID into numerical
biomarker$PatientID <- as.numeric(biomarker$PatientID)

# Create a new dataframe to join up the biomarker and covariates dataset based on PatientID
fulldata <- full_join(covariates,biomarker)

# Select column names for biomarker levels at inclusion (time = 0)
biomarker_name_inclusion <- colnames(fulldata[,7:15])

# Visualise the VAS level by biomarker level to have an preliminary data discovery
fulldata %>%
  drop_na(any_of(biomarker_name_inclusion))%>%
  select("PatientID","VAS-at-inclusion","IL-8_0weeks":"CSF-1_0weeks")%>%
  pivot_longer(names_to="BiomarkerName",values_to="BiomarkerLevel",-c("PatientID","VAS-at-inclusion"))%>%
  ggplot(aes(x=`VAS-at-inclusion`,y=BiomarkerLevel))+
  geom_point()+
  facet_wrap(~BiomarkerName)+
  ggtitle("VAS at inclusion by Biomarker Level")
         
# Define two datasets based on the VAS level at inclusion - high VAS (>=5) and low VAS (<5). If the patient has missing values for biomarker level at inclusion, remove the record accordingly.
highVASdata <- fulldata %>%
  filter(`VAS-at-inclusion`>= 5) %>%
  drop_na(any_of(biomarker_name_inclusion))
lowVASdata <- fulldata %>%
  filter(`VAS-at-inclusion`<5) %>%
  drop_na(any_of(biomarker_name_inclusion))

# Compare the descriptive statistics for the two groups
highVASstat <- highVASdata %>%
  select("IL-8_0weeks":"CSF-1_0weeks")%>%
  describe()%>%
  mutate(VASgroup="HighVAS")

highVASstat[,1] <- rownames(highVASstat)

lowVASstat <- lowVASdata %>%
  select("IL-8_0weeks":"CSF-1_0weeks")%>%
  describe()%>%
  mutate(VASgroup="LowVAS")
lowVASstat[,1] <- rownames(lowVASstat)

# Joint the descriptive statistics of high VAS and low VAS groups into one data frame for easier computation later
overallVASstat <- full_join(highVASstat,lowVASstat)%>%
  rename(nameofbiomarker=vars)%>%
  relocate(VASgroup,.after=nameofbiomarker)%>%
  select("nameofbiomarker","VASgroup","min","mean","median","max","sd")
# Visualise the differences in descriptive statistics
overallVASstat %>%
  pivot_longer(names_to="Descriptive_Statistics",values_to="Values",-c(nameofbiomarker,VASgroup))%>%
  ggplot(aes(x=Values,y=nameofbiomarker))+
  geom_line()+
  geom_point(aes(color=VASgroup))+
  theme(legend.position = "top")+
  facet_wrap(~Descriptive_Statistics,scales="free_x")+
  ggtitle("Comparison of Descriptive Statistics for High and Low VAS Groups by Biomarker")
# Convert the overallVASstat from long to wider format for easier comparison
overallVASstat_tbl <- overallVASstat %>%
  pivot_wider(names_from=VASgroup,values_from=all_of(c("min","mean","median","max","sd")))   

# Visualise the distribution of biomarker values by VAS group
highVASdata %>%
  select("PatientID","IL-8_0weeks":"CSF-1_0weeks")%>%
  pivot_longer(names_to="BiomarkerName",values_to="BiomarkerLevel",-"PatientID")%>%
  ggplot(aes(x=BiomarkerLevel))+
  geom_histogram(aes(y=after_stat(density)))+
  geom_density()+
  facet_wrap(~BiomarkerName,scales="free")+
  ggtitle("Biomarker Level Distribution for High VAS Group")

lowVASdata %>%
  select("PatientID","IL-8_0weeks":"CSF-1_0weeks")%>%
  pivot_longer(names_to="BiomarkerName",values_to="BiomarkerLevel",-"PatientID")%>%
  ggplot(aes(x=BiomarkerLevel))+
  geom_histogram(aes(y=after_stat(density)))+
  geom_density()+
  facet_wrap(~BiomarkerName,scales="free")+
  ggtitle("Biomarker Level Distribution for Low VAS Group")


# Conduct both uncorrected and corrected t-tests for each of the biomarker using for loop and save the p-value results into the dataframe t_test_result

t_test_result <- data.frame()

for (i in 7:15){
  biomarker_t_test <- t.test(highVASdata[,i],lowVASdata[,i],alternative="two.sided",var.equal=FALSE)
  t_test_result[1,i-6] <- biomarker_t_test$p.value
  t_test_result[2,i-6] <- if(biomarker_t_test$p.value > 0.05){
    "Do Not Reject"} else {
      "Reject"
    }
  t_test_result[3,i-6] <- if(biomarker_t_test$p.value > 0.05/9){
    "Do Not Reject"} else {
      "Reject"
    }
}
colnames(t_test_result) <- colnames(highVASdata[,7:15])
rownames(t_test_result) <- rbind("p-value of t-test","Pre-Bonferroni correction: Reject Null Hypothesis?", "Post-Bonferroni correction: Reject Null Hypothesis?")
                             
############ Question (2) ############

# Cleanse the data and filter all any patient with missing data for time = 0 and 12 weeks

regdata <- fulldata %>%
  select("Age":"CSF-1_0weeks")%>%
  drop_na()%>%
  relocate(`Vas-12months`,.before=Age)%>%
  relocate(`VAS-at-inclusion`,.after=`Vas-12months`)

# Set training and testing datasets
set.seed(1332)
sample_id <- sample(1:115,size = 115*0.8)
regdata_train <- regdata[sample_id,]
regdata_test <- regdata[-sample_id,]

# Train the data through regression
regmodel <- lm(`Vas-12months`~.,data=regdata_train)

# Review the results and visualise the residuals
summary(regmodel)
ggplot(regdata_train,aes(x=regmodel$residuals))+
  geom_histogram(aes(y=after_stat(count)),colour = "black",fill="light blue")+
  geom_density(aes(y=after_stat(density)*40),colour = "red",bins=20)+
  scale_y_continuous(name="count",sec.axis=sec_axis(transform = ~./40, name = "density"))+
  ggtitle("Distribution of Regression Model Residuals")+
  theme(axis.title.y.right = element_text(color = "red"))+
  theme(axis.text.y.right = element_text(color = "red"))
# Note the density curve and its axis are scaled at opposite direction in order to make the curve more visible and making sure the axis reflects its values
par(mfrow=c(1,2))
plot(predict(regmodel),residuals(regmodel),main="Residuals by Fitted Values")
qqnorm(regmodel$residuals)
qqline(regmodel$residuals,col="red")
boxplot(regmodel$residuals,main="Boxplot for Residuals",ylim=c(-6,6))


# Evaluate the model using the testing data. The testing data is fitted into the model and the predicted values and the differences with actual values are captured in the regdata_test dataset
regmodeltest <- predict(regmodel,regdata_test)
regdata_test[,15] <- regmodeltest
regdata_test[,16] <- regdata_test[,1]-regdata_test[,15]
colnames(regdata_test)[15:16] <- c("Predicted Vas-12months","Prediction Error")

# Visualise the prediction errors
regdata_test %>%
  ggplot(aes(x=`Prediction Error`))+
  geom_histogram(aes(y=after_stat(count)),colour = "black",fill="light blue",bins=10)+
  geom_density(aes(y=after_stat(density)*15),colour = "red")+
  scale_y_continuous(name="count",sec.axis=sec_axis(transform = ~./15, name = "density"))+
  ggtitle("Distribution of Prediction Error")+
  theme(axis.title.y.right = element_text(color = "red"))+
  theme(axis.text.y.right = element_text(color = "red"))
# Note the density curve and its axis are scaled at opposite direction in order to make the curve more visible and making sure the axis reflects its values
boxplot(regdata_test[,16],main="Boxplot for Prediction Error",ylim=c(-6,6))

