#Author: Ivo Pieniak
#Exploratory analysis of RNAseq data 
#Dataset: GDS5093
#PubMed : "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4116428/'



#Importing the library GEOquery
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
#Loading the library
library(GEOquery)
# Importing data from GEO
dengue_data = getGEO("GDS5093", GSEMatrix = TRUE)
dengue_table = Table(dengue_data)

gene_names= as.character(dengue_table$IDENTIFIER)
head(gene_names)

#GSE conversion into ExpressionSet object
#assayData: hints methods used to access different data components
#phenoData: meta-data describing samples 
#featureData: annotations and chip/tech features used 
#experimentData: flexible structure to describe the experiment
eset = GDS2eSet(deng_virus, do.log2= TRUE)
eset 

#Summary of the phenotypic data; table generation highlighting 4 population count
# Covalescent, dengue fever, hemorrhagic fever, control
pdat = pData(eset)
pclass = as.factor(pdat$disease.state)
names(pclass) = sampleNames(eset)
table(pclass)
levels(pclass)[levels(pclass)=="Dengue Fever"] <- "infected_dengue_fev"
levels(pclass)[levels(pclass)=="Dengue Hemorrhagic Fever"] <- "infected_dengue_hemorr_fev"
levels(pclass)[levels(pclass)=="Convalescent"] <- "convalescent"
levels(pclass)[levels(pclass)=="healthy control"] <- "control"
table(pclass)

# labelling the rows with the gene_names 
expression_set = exprs(eset)
rownames(expression_set) = gene_names
summary(expression_set)

# Setting up the indices for each of the investigated groups of samples
convalescentIndx = which ((pclass == "convalescent") == TRUE)
infected_dengue_fevIndx = which ((pclass == "infected_dengue_fev") == TRUE)
infected_dengue_hemorr_fevIndx = which ((pclass == "infected_dengue_hemorr_fev") == TRUE)
controlIndx = which ((pclass == "control") == TRUE)

colour = NULL
colour[controlIndx] = "blue"
colour[convalescentIndx] = "green"
colour[infected_dengue_fevIndx] = "red"
colour[infected_dengue_hemorr_fevIndx] = "orange"
#Presentation of the expression_set indicating different populations
boxplot(expression_set, col=colour, las = 2)

############# Principal Component Analysis ###########
# Singular Value decomposistion 

expression_pca = prcomp(t(expression_set), scale = TRUE)
attributes(expression_pca)
s= summary(expression_pca)
str(s)
# Explained Variance 
expVar = s$importance[2,] * 100
expVar 
barplot(expVar, xlab = "Principal Components", ylab = "Explained Variance (%)",
        col ="steelblue")
#Cumulative Variance
cumVar = s$importance[3,] * 100
cumVar
plot(cumVar, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance",
     xlab="Principal Components")

#Assesment of the PCA scores and loadings 
Xscores <- expression_pca$x
#plot(Xscores, xlab="PC1", ylab="PC2")
plot(Xscores, xlab="PC1", ylab="PC2", pch=21, bg=colour, cex=0.7,
     cex.lab=0.7, cex.axis = 0.7)

############# Differential Gene Expression analysis ############
###biocLite limma package
#biocLite("limma")
library(limma)

design = model.matrix(~0+pclass)
head(design)
# Creating a model
fit = lmFit( expression_set, design)
# Creating comparison pairs for makeContrasts function for further stats analysis
contrasts = makeContrasts(pclasscontrol - pclassconvalescent,
                          pclasscontrol - pclassinfected_dengue_fev,
                          pclasscontrol - pclassinfected_dengue_hemorr_fev,
                          pclassconvalescent - pclassinfected_dengue_fev,
                          pclassconvalescent - pclassinfected_dengue_hemorr_fev,
                          pclassinfected_dengue_fev - pclassinfected_dengue_hemorr_fev,
                          levels = design)

fit = contrasts.fit(fit,contrasts)
fit = eBayes(fit)
#Stats analysis table including p-values, F-values for volcano plots
toptable = topTable(fit, resort.by="p", number = 250)  
head(toptable)

hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', main='adjust P-value distribution for differentially expressed genes ', xlab="adjust p-value")
# number of differentailly expressed genes, 1 = upregulation, -1 = downregulation
table(decideTests(fit)@.Data)
#Volcano plot
with(toptable, plot(pclasscontrol...pclassconvalescent, -log10(P.Value), pch= 20, col = 'blue',
                    main = 'Volcano Plot: control&convalescent', xlim = c(-1,1)))
with(toptable, plot(pclasscontrol...pclassinfected_dengue_fev, -log10(P.Value), pch= 20, col = 'blue',
                    main = 'Volcano Plot: control&infected dengue fever', xlim = c(-1,1)))
with(toptable, plot(pclasscontrol...pclassinfected_dengue_hemorr_fev, -log10(P.Value), pch= 20, col = 'blue', 
                    main = 'Volcano Plot: control&infected dengue hemorrhagic fever', xlim = c(-1,1)))
with(toptable, plot(pclassconvalescent...pclassinfected_dengue_fev, -log10(P.Value), pch= 20, col = 'blue',
                    main = 'Volcano Plot: convalescent&infected dengue fever', xlim = c(-1,1)))
with(toptable, plot(pclassconvalescent...pclassinfected_dengue_hemorr_fev, -log10(P.Value), pch= 20, col = 'blue', 
                    main = 'Volcano Plot: convalescent& infected dengue hemorrhagic fever', xlim = c(-1,1)))
with(toptable, plot(pclassinfected_dengue_fev...pclassinfected_dengue_hemorr_fev , -log10(P.Value), pch= 20, col = 'blue',
                    main = 'Volcano Plot: dengue fever & dengue hemorrhagic fever', xlim = c(-1,1)))
