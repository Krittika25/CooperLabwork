setwd("~/Documents/TEfasta")
options(stringsAsFactors=FALSE)
library(ggplot2)
#install.packages("hrbrthemes")
library(hrbrthemes)
#library(gridExtra)
### ggpubr contains the ggarrange function
library(ggpubr)

### Read in the new TE age estimate with positions file
### that LC created with the Perl script saveTEpos.pl
leotidata = read.table("TEpos/leoti_TEpos.txt", header=TRUE)
riodata = read.table("TEpos/rio_TEpos.txt",header = TRUE)
cadata = read.table("TEpos/CA_TEpos.txt", header = TRUE)
pi229data = read.table("TEpos/pi229841_TEpos.txt", header = TRUE)
pi297data = read.table("TEpos/pi297155_TEpos.txt", header = TRUE)
pi300data = read.table("TEpos/pi300119_TEpos.txt", header = TRUE)
pi329data = read.table("TEpos/pi329311_TEpos.txt", header = TRUE)
pi506data = read.table("TEpos/pi506069_TEpos.txt", header = TRUE)
pi510data = read.table("TEpos/pi510757_TEpos.txt", header = TRUE)
pi655data = read.table("TEpos/pi655972_TEpos.txt", header = TRUE)

### Sort the table by chromosome and position
leotidata2 = leotidata[order(leotidata[,1], leotidata[,7]),]
riodata2 = riodata[order(riodata[,1], riodata[,7]),]
cadata2 = cadata[order(cadata[,1], cadata[,7]),]
pi229data2 = pi229data[order(pi229data[,1], pi229data[,7]),]
pi297data2 = pi297data[order(pi297data[,1], pi297data[,7]),]
pi300data2 = pi300data[order(pi300data[,1], pi300data[,7]),]
pi329data2 = pi329data[order(pi329data[,1], pi329data[,7]),]
pi506data2 = pi506data[order(pi506data[,1], pi506data[,7]),]
pi510data2 = pi510data[order(pi510data[,1], pi510data[,7]),]
pi655data2 = pi655data[order(pi655data[,1], pi655data[,7]),]

### Get a vector of all possible 100kb windows between position 1 and chrEnd
win.starts = seq(from=1, to=80884392, by=100000)

### Create a results table to hold the average TE age for every window on every chromosome
CHR = c(rep(1, length(win.starts)),
        rep(2, length(win.starts)),
        rep(3, length(win.starts)),
        rep(4, length(win.starts)),
        rep(5, length(win.starts)),
        rep(6, length(win.starts)),
        rep(7, length(win.starts)),
        rep(8, length(win.starts)),
        rep(9, length(win.starts)),
        rep(10, length(win.starts)))
POS = rep(win.starts, 10)
MEAN_AGE = rep(0, length(POS))

leotiresults = data.frame(CHR, POS, MEAN_AGE)
rioresults = data.frame(CHR, POS, MEAN_AGE)
caresults = data.frame(CHR, POS, MEAN_AGE)
pi229results = data.frame(CHR, POS, MEAN_AGE)
pi297results = data.frame(CHR, POS, MEAN_AGE)
pi300results = data.frame(CHR, POS, MEAN_AGE)
pi329results = data.frame(CHR, POS, MEAN_AGE)
pi506results = data.frame(CHR, POS, MEAN_AGE)
pi510results = data.frame(CHR, POS, MEAN_AGE)
pi655results = data.frame(CHR, POS, MEAN_AGE)

### For every row in the results file, get the subset of data corresponding to that chromosome and window,
### calculate and save the mean TE age
for (i in 1:nrow(leotiresults)) {
  sub.data = subset(leotidata2, leotidata2$Chromosome==leotiresults$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=leotiresults$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (leotiresults$POS[i] + 100000))
  leotiresults$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(rioresults)) {
  sub.data = subset(riodata2, riodata2$Chromosome==rioresults$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=rioresults$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (rioresults$POS[i] + 100000))
  rioresults$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(caresults)) {
  sub.data = subset(cadata2, cadata2$Chromosome==caresults$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=caresults$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (caresults$POS[i] + 100000))
  caresults$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi229results)) {
  sub.data = subset(pi229data2, pi229data2$Chromosome==pi229results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi229results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi229results$POS[i] + 100000))
  pi229results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi297results)) {
  sub.data = subset(pi297data2, pi297data2$Chromosome==pi297results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi297results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi297results$POS[i] + 100000))
  pi297results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi300results)) {
  sub.data = subset(pi300data2, pi300data2$Chromosome==pi300results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi300results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi300results$POS[i] + 100000))
  pi300results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi329results)) {
  sub.data = subset(pi329data2, pi329data2$Chromosome==pi329results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi329results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi329results$POS[i] + 100000))
  pi329results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi506results)) {
  sub.data = subset(pi506data2, pi506data2$Chromosome==pi506results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi506results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi506results$POS[i] + 100000))
  pi506results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi510results)) {
  sub.data = subset(pi510data2, pi510data2$Chromosome==pi510results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi510results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi510results$POS[i] + 100000))
  pi510results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}

for (i in 1:nrow(pi655results)) {
  sub.data = subset(pi655data2, pi655data2$Chromosome==pi655results$CHR[i])
  sub.data = subset(sub.data, sub.data$Position>=pi655results$POS[i])
  sub.data = subset(sub.data, sub.data$Position < (pi655results$POS[i] + 100000))
  pi655results$MEAN_AGE[i] = mean(sub.data$Age_MY)
}


### Replace NAs with Zeros
leotiresults$MEAN_AGE[is.na(leotiresults$MEAN_AGE)] = 0
rioresults$MEAN_AGE[is.na(rioresults$MEAN_AGE)] = 0
caresults$MEAN_AGE[is.na(caresults$MEAN_AGE)] = 0
pi229results$MEAN_AGE[is.na(pi229results$MEAN_AGE)] = 0
pi297results$MEAN_AGE[is.na(pi297results$MEAN_AGE)] = 0
pi300results$MEAN_AGE[is.na(pi300results$MEAN_AGE)] = 0
pi329results$MEAN_AGE[is.na(pi329results$MEAN_AGE)] = 0
pi506results$MEAN_AGE[is.na(pi506results$MEAN_AGE)] = 0
pi510results$MEAN_AGE[is.na(pi510results$MEAN_AGE)] = 0
pi655results$MEAN_AGE[is.na(pi655results$MEAN_AGE)] = 0


### Plot the results by chromosome (this example does chromosome 1)
leotichr.results = subset(leotiresults, leotiresults$CHR==1)
riochr.results = subset(rioresults, rioresults$CHR==1)
cachr.results = subset(caresults, caresults$CHR==1)
pi229chr.results = subset(pi229results, pi229results$CHR==1)
pi297chr.results = subset(pi297results, pi297results$CHR==1)
pi300chr.results = subset(pi300results, pi300results$CHR==1)
pi329chr.results = subset(pi329results, pi329results$CHR==1)
pi506chr.results = subset(pi506results, pi506results$CHR==1)
pi510chr.results = subset(pi510results, pi510results$CHR==1)
pi655chr.results = subset(pi655results, pi655results$CHR==1)




### Alt color scheme 1

p1 <- ggplot(data=leotichr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p2 <- ggplot(data=riochr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()


p3 <- ggplot(data=cachr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p4 <- ggplot(data=pi229chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p5 <- ggplot(data=pi297chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p6 <- ggplot(data=pi300chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p7 <- ggplot(data=pi329chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p8 <- ggplot(data=pi506chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p9 <- ggplot(data=pi510chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()

p10 <- ggplot(data=pi655chr.results[,2:3], aes(x=POS, y=1, fill=MEAN_AGE)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme_ipsum()


#grid.arrange(p1,p2,p3,p4,p5, nrow=5)

myplot <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow=10,common.legend = TRUE, labels= c("Leoti","Rio","CA","pi229841","pi297155","pi300119","pi319311","pi506069","pi510757","pi655972"), legend ="bottom")

ggsave("TEagedistplot.png",myplot, width = 8,height = 20)
