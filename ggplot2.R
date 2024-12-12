#Demo ####
#Graphics: ggplot() ####
library(ggplot2)
Dex1 <- data.frame(drug=rep(c("A","B"),6),
                   conc=c(0.7,1.6,0.4,1.7,0.5,1.4,0.8,1.8,1.5,
                          2.1,1.6,0.2))

ggplot(Dex1,aes(drug,conc))
#you just define the data for the plot
#now you have to define the look of the plot
ggplot(Dex1,aes(drug,conc)) + 
  stat_summary(fun=mean,geom="bar",width=0.2)

#2 factors
Dex1$gender <- c("F","M","F","F","M","M","F","F","F","F","F","F")

ggplot(Dex1,aes(drug,conc,fill=gender))+
  stat_summary(fun=mean,geom='bar',position='dodge')

ggplot(Dex1,aes(drug,conc,colour=gender))+
  stat_summary(fun=mean,geom='bar',position='dodge')

#color for geoms you cannot fill
ggplot(Dex1,aes(drug,conc,color=gender)) + geom_point()
#shape for geoms that can have a shape
ggplot(Dex1,aes(drug,conc,shape=gender)) + geom_point()
#shape for geoms that can have a size
ggplot(Dex1,aes(drug,gender,size=conc)) + geom_point()

#Plotting statistics ####
#Bars representing a mean
#you can save the data of a plot in a variable
p <- ggplot(Dex1,aes(drug,conc))

#Bars representing a median
p + stat_summary(fun=median,geom="bar",width=1)

p <- ggplot(Dex1,aes(drug,conc)) +
  stat_summary(fun=mean,geom="bar",width=0.3)
#Add error bars that represent 95% CI
p + stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2) +
  ylim(0,3)

#Error bars that represent sem
p + stat_summary(fun.data=mean_se,geom="errorbar",width=0.2) 

#Error bars represent 2SDs
p + stat_summary(fun.data=mean_sdl,geom="errorbar",width=0.2) 
?mean_sdl
#Error bars represent 1 SD
p + stat_summary(fun.data=mean_sdl,geom="errorbar",
                 fun.args=list(mult=1),width=0.2) 

#Error bars represent IQR
?median_hilow
#p + stat_summary(fun.data=median_hilow,geom="errorbar",
#                 fun.args=list(?=0.5),width=0.1)

#Why factors are important
ggplot(mtcars,aes(cyl,mpg)) + stat_summary(fun=mean,geom="bar")
ggplot(mtcars,aes(factor(cyl),mpg)) + stat_summary(fun=mean,geom="bar")

#Position of bars in plots with 2 factors
p <- ggplot(mtcars,aes(factor(cyl),mpg,fill=factor(am)))
p + stat_summary(fun=mean,geom="bar")
p + stat_summary(fun=mean,geom="bar",position="dodge")
p + stat_summary(fun=mean,geom="bar",position="stack")

#Stacked bar chart where the bars are not statistics
data <- read.delim("Rdata/Bacteria.txt")
#If you don't want to plot statistics: fun=identity
ggplot(data,aes(Sample,RA,fill=Phylum)) + 
  stat_summary(fun=identity,geom="bar",position="stack",width=0.5) 

#scatter plot ####
ggplot(Dex1,aes(drug,conc)) + geom_point() 
ggplot(Dex1,aes(drug,conc)) + geom_point(position="jitter") 
ggplot(Dex1,aes(drug,conc)) + geom_jitter()
?geom_jitter
ggplot(Dex1,aes(drug,conc)) + geom_jitter(width=0.1)

#beeswarm plots
library(ggbeeswarm)
babies <- read.csv2("Rdata/babies.csv")
ggplot(babies,aes(factor(smoke),wt)) + geom_point()
ggplot(babies,aes(factor(smoke),wt)) + geom_quasirandom()

#Plot with points, bars and error bars
#R plots in the order that you define layers
p + geom_point() + stat_summary(fun=mean,geom="bar") +
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.25)

#first bars then error bars then data
p + stat_summary(fun=mean,geom="bar") +
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.25) + 
  geom_point()

#Bar chart with errorbars and points for 2 factors
p <- ggplot(mtcars,aes(factor(cyl),mpg,fill=factor(am))) 
p + stat_summary(fun=mean,geom="bar",position="dodge") + 
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=0.2) +
  geom_point()

#dodge is different for every layer
p + stat_summary(fun=mean,geom="bar",position="dodge") + 
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",
               position="dodge",width=0.2) +
  geom_point(position="dodge") 

#Define a global dodge for error bars and dots
d <- position_dodge(width=0.9)
p + stat_summary(fun=mean,geom="bar",position=d,width=0.5) + 
  stat_summary(fun.data=mean_cl_normal,geom="errorbar",
               position=d,width=0.2) +
  geom_point(position=d)

#Box plot ####
ggplot(mtcars,aes(cyl,mpg)) + geom_boxplot()
ggplot(mtcars,aes(factor(cyl),mpg)) + geom_boxplot()

#Horizontal plot
ggplot(mtcars,aes(mpg,factor(cyl))) + geom_boxplot()
#Discretize X axis
ggplot(mtcars,aes(cyl,mpg)) + geom_boxplot(aes(cut_width(cyl,width=1)))

?geom_boxplot
p <- ggplot(mtcars,aes(factor(cyl),mpg)) 

#Change colors/shapes
p + geom_boxplot(outlier.color="red",outlier.shape=1)

#Make box width proportial to n
p + geom_boxplot(varwidth=TRUE)

#Make notched boxplot
p + geom_boxplot(notch=TRUE)
#https://sites.google.com/site/davidsstatistics/home/notched-box-plots

#Box plot with points
p + geom_boxplot() + geom_jitter(width=0.05)

#Add mean to boxplot
p + geom_boxplot() + 
  stat_summary(fun=mean,geom="point",color="red",size=3)

#Box plot with 2 factors
p + geom_boxplot(aes(color=factor(am)))

#Box plot with error bars
p + stat_boxplot(geom="errorbar",width=0.4) + geom_boxplot(width=0.4)

#Violin plots ####
p + geom_violin()
?geom_violin
#Set width proportional to n?
p + geom_violin(scale="count")

#Histogram ####
ggplot(mtcars,aes(mpg)) + geom_histogram()
?geom_histogram

#Adding lines to a plot ####
ggplot(mtcars,aes(wt,-mpg)) + geom_point(aes(color=factor(am))) +
  geom_abline(intercept=-35,slope=5) + xlim(0,6)
#This is not a regression line !

#Regression lines
#method=lm : linear regression -> straight lines
#y=ax+b (normal distribution)

#Adding regression lines for 2 factors
p <- ggplot(mtcars,aes(wt,mpg,color=factor(am))) + geom_point() 
p + geom_smooth(method="lm")
?geom_smooth

#method=glm : generalised linear regression -> straight lines
#y=f(ax) (no need for normal distribution - f depends on distribution 
p + geom_smooth(method="glm")
#https://www.researchgate.net/post/What_is_the_difference_between_the_general_linear_model_GLMand_generalized_linear_model_GZLM

#method=loess : nonlinear regression -> curves
#made by succession of linear regressions in small moving window
p + geom_smooth()

#method=gam : generalised additive model -> curves
#y=f(x1) + f(x2) + b
p + geom_smooth(method="gam")
#you can choose models using method + formula

#Add formula and R2 to the plot
library(ggpmisc)
myformula <- y ~ x
p + geom_smooth(method="lm") + 
  stat_poly_eq(formula=myformula,parse=TRUE,
               aes(label=paste(..eq.label..,..rr.label..,sep="~~~")))

#Text geoms####
mtcars$car <- rownames(mtcars)
mtcar <- subset(mtcars,am==1)
ggplot(mtcar,aes(wt,mpg)) + geom_point(color="red") +
  geom_text(aes(label=car))

library(ggrepel)
ggplot(mtcar,aes(wt,mpg)) + geom_point(color="red") + 
  geom_text_repel(aes(label=car))

#Scales####
#change y axis
ggplot(mtcars,aes(x=factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge") +
  scale_y_continuous(breaks=c(0,5,10,15,20,25),name="test")

#change x axis
ggplot(mtcars,aes(x=factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge") +
  scale_x_discrete(breaks=c("4","8"))

#plot with 2 Y-axes: one right, one left  ####
d <- data.frame(time=1:6,meas1=1:6,meas2=c(10,25,30,25,20,10))
#R needs a scaling between the 2 Y-axes
scaleFactor <- max(d$meas1)/max(d$meas2)
ggplot(d,aes(x=time)) +
  geom_line(aes(y=meas1),color="blue") +
  geom_line(aes(y=meas2*scaleFactor),color="red") +
  scale_y_continuous(name="conc1",
                     sec.axis=sec_axis(~./scaleFactor,name="conc2"))

#change colors
ggplot(mtcars,aes(x=factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge") + 
  scale_fill_manual(name="transmission",
                    labels=c("auto","man"),
                    values=c("red","blue"))

#change shapes
ggplot(mtcars,aes(wt,mpg,shape=factor(cyl),color=factor(cyl))) +
  geom_point(size=5) +
  scale_shape_manual(values=c("\u25C6","\u25D8","\u25BC"))
# where do you get the codes for the shapes ?
# https://jrgraphix.net/r/Unicode/25A0-25FF

#Create your own palette
library(RColorBrewer)
#Show all the color schemes available
display.brewer.all()
#Use a palette from RColorBrewer
?scale_fill_brewer
#Blues is the default palette
p <- ggplot(mtcars,aes(factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge",width=0.5) 
p + scale_fill_brewer()
#Change palette 
p + scale_fill_brewer(palette="Set3")

#Make a large color range using colorRampPalette()
#Purples palette only contains 9 colors, suppose I need 27 colors
#Start from Purples palette and use all 9 colors from the palette
purples <- brewer.pal(9,"Purples")
#Create a color range yourself
purple_range <- colorRampPalette(purples)
#Pick 27 colors from your color range for the plot
unique(mtcars$disp)
ggplot(mtcars,aes(factor(am),mpg,color=factor(disp))) + geom_point() + 
  scale_color_manual(name="disp",values=purple_range(27))

#Themes####
#Built in themes####
ggplot(mtcars,aes(factor(gear),mpg)) + 
  geom_point() + theme_bw()

#Change the background color of the legend
ggplot(mtcars,aes(factor(cyl),mpg,color=factor(am))) + 
  geom_point() +
  theme(legend.key=element_rect(fill="white",color="black"))

#part of the title italic
my_x_title <- expression(paste("No. of ",italic("cylinders")))
my_title <- paste("No. of ",italic("cylinders"))
ggplot(mtcars,aes(factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge") + 
  labs(x=my_x_title)

#Greek letters in title ####
my_x_title <- expression(paste("Symbol for alpha ",alpha))
ggplot(mtcars,aes(factor(cyl),mpg,fill=factor(am))) + 
  stat_summary(fun=mean,geom="bar",position="dodge") + 
  labs(x=my_x_title)
