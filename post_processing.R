## Summarizing USDOS output into results figures ##
## Code by Deedra Murrieta ##

#### Setup ####
library(stats); library(maps); library(ClayMapsPackage); library(mapdata); library(tidyr); library(fields); library(foreach);
library(dplyr); library(RColorBrewer); library(rgdal); library(reshape2); library(reshape); library(data.table); library(knitr);
library(ggpubr); library(ggplot2); library(ggmap); library(kableExtra); library(magrittr);library(zoo);library(gridExtra)

# Setting the working directory
path <- "C:/USDOS/"
setwd(path)

# Providing color palettes for maps and generic scales
color_red <- brewer.pal(8, "OrRd"); color_purple <- brewer.pal(8, "Purples")
color_orange <- brewer.pal(8, "Oranges"); color_bluepurple <- brewer.pal(8, "BuPu")
color_diverging <- brewer.pal(8, "RdBu"); color_blue <- brewer.pal(8, "Blues")
cbPalette <- c("#D55E00", "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#000000", "#F0E442", "#E69F00")
color_paired <- c("#000000", "#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                  "#CAB2D6", "#6A3D9A")
color_paired_ordered <- c()
color_YlOrRd <- rev(c("#000000", "#700000", "#8F0000", "#B10000", "#CE2900", "#E25700", "#F37B00", 
                      "#FF9C08", "#FFB954", "#FFD47F", "#FFEBA3", "#FFFBBF"))

scale_100 <- seq(0, 100, 10)

## Find all of the summmary files in the working directory
summary.fnames <- list.files(path = getwd(), recursive = TRUE, pattern = "_summary.txt", full.names = FALSE)
summary.fnames.90 <- summary.fnames[grep("90|75|Base|base", summary.fnames)]

## Find all the detail files in the working directory
detail.fnames <- list.files(path = path, recursive = TRUE, pattern = "_detail.txt", full.names = FALSE)
detail.fnames.90 <- detail.fnames[grep("90|75|NoCtrl|noctrl|base", detail.fnames)]

#### Duration of Infection (Summary Files) ####

county.dur <- county.fips

# Duration of infection

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "Duration"))
  summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
  summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
  summary.res <- as.data.frame(summary.res)
  county.dur <- merge(county.dur, summary.res, by.x = "fips", by.y = "Seed_FIPS", all=TRUE)
}
names(county.dur)[seq(4,ncol(county.dur),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.dur))/2-1))))
names(county.dur)[seq(3,ncol(county.dur),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.dur))/2-1))))

write.csv(county.dur, "county_dur.csv")
#county.dur <- read.csv("county_dur.csv", header = TRUE)
#county.dur <- county.dur[,-1]

dup_fips <- rownames(county.dur)[duplicated(county.dur$fips) == TRUE]
dup_name <- rownames(county.dur)[duplicated(county.dur$polyname) == TRUE]
county.dur <- county.dur[!(rownames(county.dur) %in% dup_fips),]
county.dur <- county.dur[!(rownames(county.dur) %in% dup_name),]

# Convert to long format for plotting
long.county.dur <- reshape(county.dur, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                           varying = list(c(grep("type_", colnames(county.dur))), c(grep("run_", colnames(county.dur)))))
rownames(long.county.dur) <- seq(1:nrow(long.county.dur))
long.county.dur <- long.county.dur[,-3]
long.county.dur <- long.county.dur[!is.na(long.county.dur$Value),]

# Rename the types of runs
long.county.dur$type <- as.factor(long.county.dur$type)
levels(long.county.dur$type)[levels(long.county.dur$type) == "base"] <- "Base"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_MvmtBan_90"] <- "IP Cull, 90% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_MvmtBan_75"] <- "IP Cull, 75% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_DC_MvmtBan_90"] <- "IP & DC Cull, 90% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_DC_MvmtBan_75"] <- "IP & DC Cull, 75% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_MvmtBan_90"] <- "IP Cull & DC Vax, 90% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_MvmtBan_75"] <- "IP Cull & DC Vax, 75% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_3km_MvmtBan_90"] <- "IP Cull & 3km Ring Vax, 90% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_3km_MvmtBan_75"] <- "IP Cull & 3km Ring Vax, 75% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_10km_MvmtBan_90"] <- "IP Cull & 10km Ring Vax, 90% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "IP_VAX_10km_MvmtBan_75"] <- "IP Cull & 10km Ring Vax, 75% Ban"
levels(long.county.dur$type)[levels(long.county.dur$type) == "base_ShipmentsOff"] <- "No Shipments"


head(long.county.dur)
summary(long.county.dur)

# Divide the data
long.county.dur.less13 <- long.county.dur[long.county.dur$Value > 13,]
long.county.dur.low <- long.county.dur[long.county.dur$Value <= 100,]
long.county.dur.high <- long.county.dur[long.county.dur$Value > 100,]

table(long.county.dur$Value)

jpeg("Dur_High_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.dur$Value[long.county.dur$Value > 100], breaks = 25, xlab = "Duration (Days)", main = expression("Duration" > "100 Days"))
dev.off()

jpeg("Dur_Low_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.dur$Value[long.county.dur$Value <= 100], breaks = 25, xlab = "Duration (Days)", main = expression("Duration" <= "100 Days"))
dev.off()

jpeg("Dur_Less13_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.dur$Value[long.county.dur$Value > 13], breaks = 25, xlab = "Duration (Days)", main = expression("Duration" > "13 Days"))
dev.off()

jpeg("Dur_Less13_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dur.less13, aes(x = long.county.dur.less13$type, y = long.county.dur.less13$Value, color = long.county.dur.less13$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") + 
  labs(x = NULL, y = "Duration of Infection", title = expression("Duration" > "13 Days")) + scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

jpeg("Dur_Low_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dur.low, aes(x = long.county.dur.low$type, y = long.county.dur.low$Value, color = long.county.dur.low$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Duration of Infection", title = expression("Duration" <= "100 Days"))
dev.off()

jpeg("Dur_High_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dur.high, aes(x = long.county.dur.high$type, y = long.county.dur.high$Value, color = long.county.dur.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Duration of Infection", title = expression("Duration" > "100 Days")) + ylim(c(100, 365))
dev.off()

jpeg("Dur_High_violin_max_swap.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dur.high, aes(x = long.county.dur.high$type, y = long.county.dur.high$Value, color = long.county.dur.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + coord_flip() +
  scale_color_manual(values = color_paired) + labs(x = NULL, y = "Duration of Infection", title = expression("Duration" > "100 Days")) + ylim(c(100, 365))
dev.off()

# Calculate the medians for each control type
NoCtrl_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "Base"])
NoCtrl_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "Base"])
IP90_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull, 90% Ban"])
IP90_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull, 90% Ban"])
IP75_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull, 75% Ban"])
IP75_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull, 75% Ban"])
IPDC90_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP & DC Cull, 90% Ban"])
IPDC90_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP & DC Cull, 90% Ban"])
IPDC75_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP & DC Cull, 75% Ban"])
IPDC75_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP & DC Cull, 75% Ban"])
IPVAX90_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & DC Vax, 90% Ban"])
IPVAX90_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & DC Vax, 90% Ban"])
IPVAX75_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & DC Vax, 75% Ban"])
IPVAX75_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & DC Vax, 75% Ban"])
IPVAX90_3_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & 3km Ring Vax, 90% Ban"])
IPVAX90_3_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & 3km Ring Vax, 90% Ban"])
IPVAX75_3_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & 3km Ring Vax, 75% Ban"])
IPVAX75_10_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & 10km Ring Vax, 75% Ban"])
IPVAX90_10_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & 10km Ring Vax, 90% Ban"])
IPVAX90_10_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & 10km Ring Vax, 90% Ban"])
IPVAX75_10_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "IP Cull & 10km Ring Vax, 75% Ban"])
IPVAX75_3_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "IP Cull & 3km Ring Vax, 75% Ban"])
NoShip_Dur_Low <- median(long.county.dur.low$Value[long.county.dur.low$type == "No Shipment"])
NoShip_Dur_High <- median(long.county.dur.high$Value[long.county.dur.high$type == "No Shipment"])


Dur_Medians_90_Low <- rbind(NoCtrl_Dur_Low, IP90_Dur_Low, IPDC90_Dur_Low, IPVAX90_Dur_Low, IPVAX90_3_Dur_Low, IPVAX90_10_Dur_Low, NoShip_Dur_Low)
Dur_Medians_90_High <- rbind(NoCtrl_Dur_High, IP90_Dur_High, IPDC90_Dur_High, IPVAX90_Dur_High, IPVAX90_3_Dur_High, IPVAX90_10_Dur_High, NoShip_Dur_High)
Dur_Medians_75_Low <- rbind(NoCtrl_Dur_Low, IP75_Dur_Low, IPDC75_Dur_Low, IPVAX75_Dur_Low, IPVAX75_3_Dur_Low, IPVAX75_10_Dur_Low, NoShip_Dur_Low)
Dur_Medians_75_High <- rbind(NoCtrl_Dur_High, IP75_Dur_High, IPDC75_Dur_High, IPVAX75_Dur_High, IPVAX75_3_Dur_High, IPVAX75_10_Dur_High, NoShip_Dur_High)

## Create individual data frames for each run type

df <- county.dur[,c(1,2,3:202)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,203:402)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,403:602)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,603:802)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,803:1002)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,1003:1202)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,1203:1402)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,1403:1602)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,1603:1802)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,1803:2002)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,2003:2202)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,2203:2402)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

df <- county.dur[,c(1,2,2403:2602)]
df1 <- df[df$type_1201 != 0 & !is.na(df$type_1201),]
newname <- df1$type_1201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_dur"), df.new)

## calculate median of each run
base_dur$median <- apply(base_dur[grep("run_", colnames(base_dur))], 1, median, na.rm = TRUE)
IP_MvmtBan_90_dur$median <- apply(IP_MvmtBan_90_dur[grep("run_", colnames(IP_MvmtBan_90_dur))], 1, median, na.rm = TRUE)
IP_MvmtBan_75_dur$median <- apply(IP_MvmtBan_75_dur[grep("run_", colnames(IP_MvmtBan_75_dur))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_90_dur$median <- apply(IP_DC_MvmtBan_90_dur[grep("run_", colnames(IP_DC_MvmtBan_90_dur))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_75_dur$median <- apply(IP_DC_MvmtBan_75_dur[grep("run_", colnames(IP_DC_MvmtBan_75_dur))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_90_dur$median <- apply(IP_VAX_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_MvmtBan_90_dur))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_75_dur$median <- apply(IP_VAX_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_MvmtBan_75_dur))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_dur$median <- apply(IP_VAX_3km_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_dur))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_dur$median <- apply(IP_VAX_3km_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_dur))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_dur$median <- apply(IP_VAX_10km_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_dur))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_dur$median <- apply(IP_VAX_10km_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_dur))], 1, median, na.rm = TRUE)
base_ShipmentsOff_dur$median <- apply(base_ShipmentsOff_dur[grep("run_", colnames(base_ShipmentsOff_dur))], 1, median, na.rm = TRUE)

# Create a data.frame with just the fips code and median to map
base_dur_med <- base_dur[,c(1,103)]
IP_MvmtBan_90_dur_med <- IP_MvmtBan_90_dur[,c(1,103)]
IP_MvmtBan_75_dur_med <- IP_MvmtBan_75_dur[,c(1,103)]
IP_DC_MvmtBan_90_dur_med <- IP_DC_MvmtBan_90_dur[,c(1,103)]
IP_DC_MvmtBan_75_dur_med <- IP_DC_MvmtBan_75_dur[,c(1,103)]
IP_VAX_MvmtBan_90_dur_med <- IP_VAX_MvmtBan_90_dur[,c(1,103)]
IP_VAX_MvmtBan_75_dur_med <- IP_VAX_MvmtBan_75_dur[,c(1,103)]
IP_VAX_3km_MvmtBan_90_dur_med <- IP_VAX_3km_MvmtBan_90_dur[,c(1,103)]
IP_VAX_3km_MvmtBan_75_dur_med <- IP_VAX_3km_MvmtBan_75_dur[,c(1,103)]
IP_VAX_10km_MvmtBan_90_dur_med <- IP_VAX_10km_MvmtBan_90_dur[,c(1,103)]
IP_VAX_10km_MvmtBan_75_dur_med <- IP_VAX_10km_MvmtBan_75_dur[,c(1,103)]
base_ShipmentsOff_dur_med <- base_ShipmentsOff_dur[,c(1,103)]


## Calculate the Upper 2.5% for each county
base_dur$upper <- apply(base_dur[grep("run_", colnames(base_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_90_dur$upper <- apply(IP_MvmtBan_90_dur[grep("run_", colnames(IP_MvmtBan_90_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_75_dur$upper <- apply(IP_MvmtBan_75_dur[grep("run_", colnames(IP_MvmtBan_75_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_90_dur$upper <- apply(IP_DC_MvmtBan_90_dur[grep("run_", colnames(IP_DC_MvmtBan_90_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_75_dur$upper <- apply(IP_DC_MvmtBan_75_dur[grep("run_", colnames(IP_DC_MvmtBan_75_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_90_dur$upper <- apply(IP_VAX_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_MvmtBan_90_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_75_dur$upper <- apply(IP_VAX_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_MvmtBan_75_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_dur$upper <- apply(IP_VAX_3km_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_dur$upper <- apply(IP_VAX_3km_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_dur$upper <- apply(IP_VAX_10km_MvmtBan_90_dur[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_dur$upper <- apply(IP_VAX_10km_MvmtBan_75_dur[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_dur))], 1, quantile, probs=0.975, na.rm = TRUE)
base_ShipmentsOff_dur$upper <- apply(base_ShipmentsOff_dur[grep("run_", colnames(base_ShipmentsOff_dur))], 1, quantile, probs=0.975, na.rm = TRUE)

# Create a data.frame with just the fips code and upper to map
base_dur_up <- base_dur[,c(1,104)]
IP_MvmtBan_90_dur_up <- IP_MvmtBan_90_dur[,c(1,104)]
IP_MvmtBan_75_dur_up <- IP_MvmtBan_75_dur[,c(1,104)]
IP_DC_MvmtBan_90_dur_up <- IP_DC_MvmtBan_90_dur[,c(1,104)]
IP_DC_MvmtBan_75_dur_up <- IP_DC_MvmtBan_75_dur[,c(1,104)]
IP_VAX_MvmtBan_90_dur_up <- IP_VAX_MvmtBan_90_dur[,c(1,104)]
IP_VAX_MvmtBan_75_dur_up <- IP_VAX_MvmtBan_75_dur[,c(1,104)]
IP_VAX_3km_MvmtBan_90_dur_up <- IP_VAX_3km_MvmtBan_90_dur[,c(1,104)]
IP_VAX_3km_MvmtBan_75_dur_up <- IP_VAX_3km_MvmtBan_75_dur[,c(1,104)]
IP_VAX_10km_MvmtBan_90_dur_up <- IP_VAX_10km_MvmtBan_90_dur[,c(1,104)]
IP_VAX_10km_MvmtBan_75_dur_up <- IP_VAX_10km_MvmtBan_75_dur[,c(1,104)]
base_ShipmentsOff_dur_up <- base_ShipmentsOff_dur[,c(1,104)]

# Create a vector for scale
dur_upper_values <- c(base_dur$upper, IP_MvmtBan_90_dur$upper, IP_MvmtBan_75_dur$upper, IP_DC_MvmtBan_90_dur$upper, 
                      IP_DC_MvmtBan_75_dur$upper, IP_VAX_MvmtBan_90_dur$upper, IP_VAX_MvmtBan_75_dur$upper, 
                      IP_VAX_3km_MvmtBan_90_dur$upper, IP_VAX_3km_MvmtBan_75_dur$upper, IP_VAX_10km_MvmtBan_90_dur$upper, 
                      IP_VAX_10km_MvmtBan_75_dur$upper, base_ShipmentsOff_dur$upper)
dur_scale <- round(quantile(dur_upper_values, probs = c(0, .21, .25, .3, .4, .5, .6, .7, .8, .9, 0.95, 1), na.rm = TRUE), 0)
dur_scale <- c(13, 16, 18, 20, 22, 24, 30, 40, 160, 285, 330, 365)

# Median maps
jpeg("Dur_Med_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Med_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_dur_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

# Upper 2.5% maps
jpeg("Dur_Up_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Dur_Up_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_dur_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = dur_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()


#### Number of Premises Infected (Summary Files) ####

# # Retrieve data from all runs
county.inf <- county.fips

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "Num_Inf"))
  summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
  summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
  summary.res <- as.data.frame(summary.res)
  county.inf <- merge(county.inf, summary.res, by.x = "fips", by.y = "Seed_FIPS", all=TRUE)
}
names(county.inf)[seq(4,ncol(county.inf),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.inf))/2-1))))
names(county.inf)[seq(3,ncol(county.inf),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.inf))/2-1))))
county.inf[nrow(county.inf)+1,] <- c(46102, "south dakota, lakota", county.inf[,3:length(county.inf)][county.inf$fips == 46007,])

dup_fips <- rownames(county.inf)[duplicated(county.inf$fips) == TRUE]
dup_name <- rownames(county.inf)[duplicated(county.inf$polyname) == TRUE]
county.inf <- county.inf[!(rownames(county.inf) %in% dup_fips),]
county.inf <- county.inf[!(rownames(county.inf) %in% dup_name),]


write.csv(county.inf, "county_inf.csv")
#county.inf <- read.csv("county_inf.csv", header = TRUE)
#county.inf <- county.inf[,-1]

# Convert to long format for plotting
long.county.inf <- reshape(county.inf, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                           varying = list(c(grep("type_", colnames(county.inf))), c(grep("run_", colnames(county.inf)))))
rownames(long.county.inf) <- seq(1:nrow(long.county.inf))
long.county.inf <- long.county.inf[,-3]
long.county.inf <- long.county.inf[!is.na(long.county.inf$Value),]

# Rename the types of runs
long.county.inf$type <- as.factor(long.county.inf$type)
levels(long.county.inf$type)[levels(long.county.inf$type) == "base"] <- "Base"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_MvmtBan_90"] <- "IP Cull, 90% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_MvmtBan_75"] <- "IP Cull, 75% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_DC_MvmtBan_90"] <- "IP & DC Cull, 90% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_DC_MvmtBan_75"] <- "IP & DC Cull, 75% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_MvmtBan_90"] <- "IP Cull & DC Vax, 90% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_MvmtBan_75"] <- "IP Cull & DC Vax, 75% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_3km_MvmtBan_90"] <- "IP Cull & 3km Ring Vax, 90% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_3km_MvmtBan_75"] <- "IP Cull & 3km Ring Vax, 75% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_10km_MvmtBan_90"] <- "IP Cull & 10km Ring Vax, 90% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "IP_VAX_10km_MvmtBan_75"] <- "IP Cull & 10km Ring Vax, 75% Ban"
levels(long.county.inf$type)[levels(long.county.inf$type) == "base_ShipmentsOff"] <- "No Shipments"

head(long.county.inf)

color_paired_NI <- c("#000000", "#B2DF8A", "#33A02C", "#A6CEE3", "#1F78B4", "#FB9A99", "#FDBF6F", "#CAB2D6", 
                     "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928" )
 
# Divide the data

long.county.inf.less10 <- long.county.inf[long.county.inf$Value > 10,]
long.county.inf.low <- long.county.inf[long.county.inf$Value <= 5000,]
long.county.inf.high <- long.county.inf[long.county.inf$Value > 5000,]

table(long.county.inf$Value)
jpeg("NI_High_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.inf$Value[long.county.inf$Value > 5000], breaks = 100, 
     main = expression("# Infected Premises" > "5000 Subset"), xlab = "# Infected Premises")
dev.off()

jpeg("NI_Low_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.inf$Value[long.county.inf$Value <= 5000 & long.county.inf$Value > 10], breaks = 100, 
     main = expression("# Infected Premises" <= "5000 Subset"), xlab = "# Infected Premises")
dev.off()

jpeg("NI_Less10_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.inf.less10, aes(x = long.county.inf.less10$type, y = long.county.inf.less10$Value, color = long.county.inf.less10$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Number of IPs",  title = expression("# Infected Premises" > "10 Subset")) + scale_color_manual(values = color_paired_NI) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

jpeg("NI_Low_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.inf.low, aes(x = long.county.inf.low$type, y = long.county.inf.low$Value, color = long.county.inf.low$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Number of IPs",  title = expression("# Infected Premises" <= "5000 Subset")) + scale_color_manual(values = color_paired_NI) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

jpeg("NI_High_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.inf.high, aes(x = long.county.inf.high$type, y = long.county.inf.high$Value, color = long.county.inf.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Number of IPs",  title = expression("# Infected Premises" > "5000 Subset")) + scale_color_manual(values = color_paired_NI) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylim(c(5000, 32000))
dev.off()

jpeg("NI_High_violin_max_swap.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.inf.high, aes(x = long.county.inf.high$type, y = long.county.inf.high$Value, color = long.county.inf.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = "Number of IPs",  title = expression("# Infected Premises" > "5000 Subset")) + scale_color_manual(values = color_paired_NI) +  
  ylim(c(5000, 32000)) + coord_flip()
dev.off()



# Calculate the medians for each control type
NoCtrl_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "Base"])
NoCtrl_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "Base"])
IP90_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull, 90% Ban"])
IP90_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull, 90% Ban"])
IP75_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull, 75% Ban"])
IP75_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull, 75% Ban"])
IPDC90_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP & DC Cull, 90% Ban"])
IPDC90_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP & DC Cull, 90% Ban"])
IPDC75_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP & DC Cull, 75% Ban"])
IPDC75_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP & DC Cull, 75% Ban"])
IPVAX90_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & DC Vax, 90% Ban"])
IPVAX90_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & DC Vax, 90% Ban"])
IPVAX75_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & DC Vax, 75% Ban"])
IPVAX75_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & DC Vax, 75% Ban"])
IPVAX90_3_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & 3km Ring Vax, 90% Ban"])
IPVAX90_3_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & 3km Ring Vax, 90% Ban"])
IPVAX75_3_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & 3km Ring Vax, 75% Ban"])
IPVAX75_10_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & 10km Ring Vax, 75% Ban"])
IPVAX90_10_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & 10km Ring Vax, 90% Ban"])
IPVAX90_10_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & 10km Ring Vax, 90% Ban"])
IPVAX75_10_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "IP Cull & 10km Ring Vax, 75% Ban"])
IPVAX75_3_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "IP Cull & 3km Ring Vax, 75% Ban"])
NoShip_NI_Low <- median(long.county.inf.low$Value[long.county.inf.low$type == "No Shipment"])
NoShip_NI_High <- median(long.county.inf.high$Value[long.county.inf.high$type == "No Shipment"])


NI_Medians_90_Low <- rbind(NoCtrl_NI_Low, IP90_NI_Low, IPDC90_NI_Low, IPVAX90_NI_Low, IPVAX90_3_NI_Low, IPVAX90_10_NI_Low, NoShip_NI_Low)
NI_Medians_90_High <- rbind(NoCtrl_NI_High, IP90_NI_High, IPDC90_NI_High, IPVAX90_NI_High, IPVAX90_3_NI_High, IPVAX90_10_NI_High, NoShip_NI_High)
NI_Medians_75_Low <- rbind(NoCtrl_NI_Low, IP75_NI_Low, IPDC75_NI_Low, IPVAX75_NI_Low, IPVAX75_3_NI_Low, IPVAX75_10_NI_Low, NoShip_NI_Low)
NI_Medians_75_High <- rbind(NoCtrl_NI_High, IP75_NI_High, IPDC75_NI_High, IPVAX75_NI_High, IPVAX75_3_NI_High, IPVAX75_10_NI_High, NoShip_NI_High)

# Create individual data frames for each run type
df <- county.inf[,c(1,2,3:202)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,203:402)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,403:602)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,603:802)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,803:1002)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,1003:1202)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,1203:1402)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,1403:1602)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,1603:1802)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,1803:2002)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,2003:2202)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,2203:2402)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)

df <- county.inf[,c(1,2,2403:2602)]
df1 <- df[df$type_1201 != 0 & !is.na(df$type_1201),]
newname <- df1$type_1201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_NI"), df.new)


## calculate median of each run
base_NI$median <- apply(base_NI[grep("run_", colnames(base_NI))], 1, median, na.rm = TRUE)
IP_MvmtBan_90_NI$median <- apply(IP_MvmtBan_90_NI[grep("run_", colnames(IP_MvmtBan_90_NI))], 1, median, na.rm = TRUE)
IP_MvmtBan_75_NI$median <- apply(IP_MvmtBan_75_NI[grep("run_", colnames(IP_MvmtBan_75_NI))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_90_NI$median <- apply(IP_DC_MvmtBan_90_NI[grep("run_", colnames(IP_DC_MvmtBan_90_NI))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_75_NI$median <- apply(IP_DC_MvmtBan_75_NI[grep("run_", colnames(IP_DC_MvmtBan_75_NI))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_90_NI$median <- apply(IP_VAX_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_MvmtBan_90_NI))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_75_NI$median <- apply(IP_VAX_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_MvmtBan_75_NI))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_NI$median <- apply(IP_VAX_3km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_NI))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_NI$median <- apply(IP_VAX_3km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_NI))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_NI$median <- apply(IP_VAX_10km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_NI))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_NI$median <- apply(IP_VAX_10km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_NI))], 1, median, na.rm = TRUE)
base_ShipmentsOff_NI$median <- apply(base_ShipmentsOff_NI[grep("run_", colnames(base_ShipmentsOff_NI))], 1, median, na.rm = TRUE)

# Create a df with just the fips code and median to map
base_NI_med <- base_NI[,c(1,103)]
IP_MvmtBan_90_NI_med <- IP_MvmtBan_90_NI[,c(1,103)]
IP_MvmtBan_75_NI_med <- IP_MvmtBan_75_NI[,c(1,103)]
IP_DC_MvmtBan_90_NI_med <- IP_DC_MvmtBan_90_NI[,c(1,103)]
IP_DC_MvmtBan_75_NI_med <- IP_DC_MvmtBan_75_NI[,c(1,103)]
IP_VAX_MvmtBan_90_NI_med <- IP_VAX_MvmtBan_90_NI[,c(1,103)]
IP_VAX_MvmtBan_75_NI_med <- IP_VAX_MvmtBan_75_NI[,c(1,103)]
IP_VAX_3km_MvmtBan_90_NI_med <- IP_VAX_3km_MvmtBan_90_NI[,c(1,103)]
IP_VAX_3km_MvmtBan_75_NI_med <- IP_VAX_3km_MvmtBan_75_NI[,c(1,103)]
IP_VAX_10km_MvmtBan_90_NI_med <- IP_VAX_10km_MvmtBan_90_NI[,c(1,103)]
IP_VAX_10km_MvmtBan_75_NI_med <- IP_VAX_10km_MvmtBan_75_NI[,c(1,103)]
base_ShipmentsOff_NI_med <- base_ShipmentsOff_NI[,c(1,103)]


## Calculate the Upper 2.5% for each county
base_NI$upper <- apply(base_NI[grep("run_", colnames(base_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_90_NI$upper <- apply(IP_MvmtBan_90_NI[grep("run_", colnames(IP_MvmtBan_90_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_75_NI$upper <- apply(IP_MvmtBan_75_NI[grep("run_", colnames(IP_MvmtBan_75_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_90_NI$upper <- apply(IP_DC_MvmtBan_90_NI[grep("run_", colnames(IP_DC_MvmtBan_90_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_75_NI$upper <- apply(IP_DC_MvmtBan_75_NI[grep("run_", colnames(IP_DC_MvmtBan_75_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_90_NI$upper <- apply(IP_VAX_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_MvmtBan_90_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_75_NI$upper <- apply(IP_VAX_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_MvmtBan_75_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_NI$upper <- apply(IP_VAX_3km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_NI$upper <- apply(IP_VAX_3km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_NI$upper <- apply(IP_VAX_10km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_NI$upper <- apply(IP_VAX_10km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_NI))], 1, quantile, probs=0.975, na.rm = TRUE)
base_ShipmentsOff_NI$upper <- apply(base_ShipmentsOff_NI[grep("run_", colnames(base_ShipmentsOff_NI))], 1, quantile, probs=0.975, na.rm = TRUE)

# Create a data.frame with just the fips code and upper to map
base_NI_up <- base_NI[,c(1,104)]
IP_MvmtBan_90_NI_up <- IP_MvmtBan_90_NI[,c(1,104)]
IP_MvmtBan_75_NI_up <- IP_MvmtBan_75_NI[,c(1,104)]
IP_DC_MvmtBan_90_NI_up <- IP_DC_MvmtBan_90_NI[,c(1,104)]
IP_DC_MvmtBan_75_NI_up <- IP_DC_MvmtBan_75_NI[,c(1,104)]
IP_VAX_MvmtBan_90_NI_up <- IP_VAX_MvmtBan_90_NI[,c(1,104)]
IP_VAX_MvmtBan_75_NI_up <- IP_VAX_MvmtBan_75_NI[,c(1,104)]
IP_VAX_3km_MvmtBan_90_NI_up <- IP_VAX_3km_MvmtBan_90_NI[,c(1,104)]
IP_VAX_3km_MvmtBan_75_NI_up <- IP_VAX_3km_MvmtBan_75_NI[,c(1,104)]
IP_VAX_10km_MvmtBan_90_NI_up <- IP_VAX_10km_MvmtBan_90_NI[,c(1,104)]
IP_VAX_10km_MvmtBan_75_NI_up <- IP_VAX_10km_MvmtBan_75_NI[,c(1,104)]
base_ShipmentsOff_NI_up <- base_ShipmentsOff_NI[,c(1,104)]

# Create a vector for scale
NI_scale <- c(base_NI$upper, IP_MvmtBan_90_NI$upper, IP_MvmtBan_75_NI$upper, IP_DC_MvmtBan_90_NI$upper, 
                      IP_DC_MvmtBan_75_NI$upper, IP_VAX_MvmtBan_90_NI$upper, IP_VAX_MvmtBan_75_NI$upper, 
                      IP_VAX_3km_MvmtBan_90_NI$upper, IP_VAX_3km_MvmtBan_75_NI$upper, IP_VAX_10km_MvmtBan_90_NI$upper, 
                      IP_VAX_10km_MvmtBan_75_NI$upper, base_ShipmentsOff_NI$upper)
NI_scale <- round(quantile(NI_scale, probs = c(0, .21, .25, .3, .4, .5, .6, .7, .8, 0.85, .9, 0.95, 0.99, 1), na.rm = TRUE), 0)
NI_scale <- c(0, 1, 2, 4, 10, 5000, 10000, 15000, 20000, 25000, 27500, 30000)

# Median maps
jpeg("NI_Med_Base_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IP90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IP75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPDC90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPDC75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX90_3_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX75_3_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX90_10_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_IPVAX75_10_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Med_ShipmentsOff_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_NI_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

# Upper 2.5% maps
jpeg("NI_Up_Base_max_Random.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IP90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IP75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPDC90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPDC75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX90_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX75_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX90_3_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX75_3_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX90_10_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_IPVAX75_10_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("NI_Up_ShipmentsOff_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_NI_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = NI_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

#### Proportion of runs with >10 infections (Summary Files) ####

## Calculate the number of >10 infections for each county
base_NI$gr10 <- apply(base_NI[grep("run_", colnames(base_NI))], 1, function(x) length(x[x>10]))
IP_MvmtBan_90_NI$gr10 <- apply(IP_MvmtBan_90_NI[grep("run_", colnames(IP_MvmtBan_90_NI))], 1, function(x) length(x[x>10]))
IP_MvmtBan_75_NI$gr10 <- apply(IP_MvmtBan_75_NI[grep("run_", colnames(IP_MvmtBan_75_NI))], 1, function(x) length(x[x>10]))
IP_DC_MvmtBan_90_NI$gr10 <- apply(IP_DC_MvmtBan_90_NI[grep("run_", colnames(IP_DC_MvmtBan_90_NI))], 1, function(x) length(x[x>10]))
IP_DC_MvmtBan_75_NI$gr10 <- apply(IP_DC_MvmtBan_75_NI[grep("run_", colnames(IP_DC_MvmtBan_75_NI))], 1, function(x) length(x[x>10]))
IP_VAX_MvmtBan_90_NI$gr10 <- apply(IP_VAX_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_MvmtBan_90_NI))], 1, function(x) length(x[x>10]))
IP_VAX_MvmtBan_75_NI$gr10 <- apply(IP_VAX_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_MvmtBan_75_NI))], 1, function(x) length(x[x>10]))
IP_VAX_3km_MvmtBan_90_NI$gr10 <- apply(IP_VAX_3km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_NI))], 1, function(x) length(x[x>10]))
IP_VAX_3km_MvmtBan_75_NI$gr10 <- apply(IP_VAX_3km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_NI))], 1, function(x) length(x[x>10]))
IP_VAX_10km_MvmtBan_90_NI$gr10 <- apply(IP_VAX_10km_MvmtBan_90_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_NI))], 1, function(x) length(x[x>10]))
IP_VAX_10km_MvmtBan_75_NI$gr10 <- apply(IP_VAX_10km_MvmtBan_75_NI[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_NI))], 1, function(x) length(x[x>10]))
base_ShipmentsOff_NI$gr10 <- apply(base_ShipmentsOff_NI[grep("run_", colnames(base_ShipmentsOff_NI))], 1, function(x) length(x[x>10]))

# For Summary Table - find the total number of times it was >10 across the counties, then divide by total possible (304900) and * 100 for percentage
Base_Prop <- sum(base_NI$gr10, na.rm = TRUE)/3049
IP90_Prop <- sum(IP_MvmtBan_90_NI$gr10, na.rm = TRUE)/3049
IPDC90_Prop <- sum(IP_DC_MvmtBan_90_NI$gr10, na.rm = TRUE)/3049
IPVAX90_Prop <- sum(IP_VAX_MvmtBan_90_NI$gr10, na.rm = TRUE)/3049
IPVAX90_3_Prop <- sum(IP_VAX_3km_MvmtBan_90_NI$gr10, na.rm = TRUE)/3049
IPVAX90_10_Prop <- sum(IP_VAX_10km_MvmtBan_90_NI$gr10, na.rm = TRUE)/3049
IP75_Prop <- sum(IP_MvmtBan_75_NI$gr10, na.rm = TRUE)/3049
IPDC75_Prop <- sum(IP_DC_MvmtBan_75_NI$gr10, na.rm = TRUE)/3049
IPVAX75_Prop <- sum(IP_VAX_MvmtBan_75_NI$gr10, na.rm = TRUE)/3049
IPVAX75_3_Prop <- sum(IP_VAX_3km_MvmtBan_75_NI$gr10, na.rm = TRUE)/3049
IPVAX75_10_Prop <- sum(IP_VAX_10km_MvmtBan_75_NI$gr10, na.rm = TRUE)/3049
NoShip_Prop <- sum(base_ShipmentsOff_NI$gr10, na.rm = TRUE)/3049

Prop_90_min <- rbind(Base_Prop, NoShip_Prop, IP90_Prop, IPDC90_Prop, IPVAX90_Prop, IPVAX90_3_Prop, IPVAX90_10_Prop)
Prop_75_min <- rbind(Base_Prop, NoShip_Prop, IP75_Prop, IPDC75_Prop, IPVAX75_Prop, IPVAX75_3_Prop, IPVAX75_10_Prop)

#### Epidemic Extent (Summary Files) ####

# Epidemic Extent

county.uni <- county.fips

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "nAffCounties"))
  summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
  summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
  summary.res <- as.data.frame(summary.res)
  county.uni <- merge(county.uni, summary.res, by.x="fips", by.y = "Seed_FIPS", all=TRUE)
}
names(county.uni)[seq(4,ncol(county.uni),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.uni))/2-1))))
names(county.uni)[seq(3,ncol(county.uni),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.uni))/2-1))))
#county.uni[nrow(county.uni)+1,] <- c(46102, "south dakota, lakota", county.uni[,3:length(county.uni)][county.uni$fips == 46007,])
dup_fips <- rownames(county.uni)[duplicated(county.uni$fips) == TRUE]
dup_name <- rownames(county.uni)[duplicated(county.uni$polyname) == TRUE]
county.uni <- county.uni[!(rownames(county.uni) %in% dup_fips),]
county.uni <- county.uni[!(rownames(county.uni) %in% dup_name),]


write.csv(county.uni, "county_uni.csv")

# Convert to long format for plotting
long.county.uni <- reshape(county.uni, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                           varying = list(c(grep("type_", colnames(county.uni))), c(grep("run_", colnames(county.uni)))))
rownames(long.county.uni) <- seq(1:nrow(long.county.uni))
long.county.uni <- long.county.uni[,-3]
long.county.uni <- long.county.uni[!is.na(long.county.uni$Value),]

# Rename the types of runs
long.county.uni$type <- as.factor(long.county.uni$type)
levels(long.county.uni$type)[levels(long.county.uni$type) == "base"] <- "Base"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_MvmtBan_90"] <- "IP Cull, 90% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_MvmtBan_75"] <- "IP Cull, 75% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_DC_MvmtBan_90"] <- "IP & DC Cull, 90% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_DC_MvmtBan_75"] <- "IP & DC Cull, 75% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_MvmtBan_90"] <- "IP Cull & DC Vax, 90% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_MvmtBan_75"] <- "IP Cull & DC Vax, 75% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_3km_MvmtBan_90"] <- "IP Cull & 3km Ring Vax, 90% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_3km_MvmtBan_75"] <- "IP Cull & 3km Ring Vax, 75% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_10km_MvmtBan_90"] <- "IP Cull & 10km Ring Vax, 90% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "IP_VAX_10km_MvmtBan_75"] <- "IP Cull & 10km Ring Vax, 75% Ban"
levels(long.county.uni$type)[levels(long.county.uni$type) == "base_ShipmentsOff"] <- "No Shipments"


head(long.county.uni)
summary(long.county.uni)

# Divide the data
long.county.uni.less10 <- long.county.uni[long.county.uni$Value > 10,]
long.county.uni.low <- long.county.uni[long.county.uni$Value <= 500,]
long.county.uni.high <- long.county.uni[long.county.uni$Value > 500,]
table(long.county.uni$Value)
hist(long.county.uni$Value[long.county.uni$Value > 2])

jpeg("Epi_High_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.uni$Value[long.county.uni$Value > 500], breaks = 25, xlab = "Epidemic Extent", main = expression("# Infected Counties" > "500 Subset"))
dev.off()

jpeg("Epi_Low_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.uni$Value[long.county.uni$Value <= 500], breaks = 25, xlab = "Epidemic Extent", main = expression("# Infected Counties" <= "500 Subset"))
dev.off()

jpeg("Epi_Less10_hist_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.uni$Value[long.county.uni$Value > 10], breaks = 25, xlab = "Epidemic Extent", main = expression("# Infected Counties" > "10 Subset"))
dev.off()

jpeg("Epi_Less10_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.uni.less10, aes(x = long.county.uni.less10$type, y = long.county.uni.less10$Value, color = long.county.uni.less10$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") + 
  labs(x = NULL, y = "Epidemic Extent", title = expression("# Infected Counties" > "10 Subset")) + scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

jpeg("Epi_Low_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.uni.low, aes(x = long.county.uni.low$type, y = long.county.uni.low$Value, color = long.county.uni.low$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Epidemic Extent", title = expression("# Infected Counties" <= "500 Subset"))
dev.off()

jpeg("Epi_High_violin_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.uni.high, aes(x = long.county.uni.high$type, y = long.county.uni.high$Value, color = long.county.uni.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) + ylim(c(500, 1600)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Epidemic Extent", title = expression("# Infected Counties" > "500 Subset"))
dev.off()

jpeg("Epi_High_violin_max_swap.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.uni.high, aes(x = long.county.uni.high$type, y = long.county.uni.high$Value, color = long.county.uni.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) +  coord_flip() + labs(x = NULL, y = "Epidemic Extent", title = expression("# Infected Counties" > "500 Subset")) + ylim(c(500, 1600))
dev.off()


## Create individual data frames for each run type
# Create individual data frames for each run type
df <- county.uni[,c(1,2,3:202)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,203:402)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,403:602)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,603:802)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,803:1002)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,1003:1202)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,1203:1402)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,1403:1602)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,1603:1802)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,1803:2002)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,2003:2202)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,2203:2402)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)

df <- county.uni[,c(1,2,2403:2602)]
df1 <- df[df$type_1201 != 0 & !is.na(df$type_1201),]
newname <- df1$type_1201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Epi"), df.new)


## calculate median of each run
base_Epi$median <- apply(base_Epi[grep("run_", colnames(base_Epi))], 1, median, na.rm = TRUE)
IP_MvmtBan_90_Epi$median <- apply(IP_MvmtBan_90_Epi[grep("run_", colnames(IP_MvmtBan_90_Epi))], 1, median, na.rm = TRUE)
IP_MvmtBan_75_Epi$median <- apply(IP_MvmtBan_75_Epi[grep("run_", colnames(IP_MvmtBan_75_Epi))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_90_Epi$median <- apply(IP_DC_MvmtBan_90_Epi[grep("run_", colnames(IP_DC_MvmtBan_90_Epi))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_75_Epi$median <- apply(IP_DC_MvmtBan_75_Epi[grep("run_", colnames(IP_DC_MvmtBan_75_Epi))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_90_Epi$median <- apply(IP_VAX_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_MvmtBan_90_Epi))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_75_Epi$median <- apply(IP_VAX_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_MvmtBan_75_Epi))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_Epi$median <- apply(IP_VAX_3km_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_Epi))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_Epi$median <- apply(IP_VAX_3km_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_Epi))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_Epi$median <- apply(IP_VAX_10km_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_Epi))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_Epi$median <- apply(IP_VAX_10km_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_Epi))], 1, median, na.rm = TRUE)
base_ShipmentsOff_Epi$median <- apply(base_ShipmentsOff_Epi[grep("run_", colnames(base_ShipmentsOff_Epi))], 1, median, na.rm = TRUE)

# Create a data.frame with just the fips code and median to map
base_Epi_med <- base_Epi[,c(1,103)]
IP_MvmtBan_90_Epi_med <- IP_MvmtBan_90_Epi[,c(1,103)]
IP_MvmtBan_75_Epi_med <- IP_MvmtBan_75_Epi[,c(1,103)]
IP_DC_MvmtBan_90_Epi_med <- IP_DC_MvmtBan_90_Epi[,c(1,103)]
IP_DC_MvmtBan_75_Epi_med <- IP_DC_MvmtBan_75_Epi[,c(1,103)]
IP_VAX_MvmtBan_90_Epi_med <- IP_VAX_MvmtBan_90_Epi[,c(1,103)]
IP_VAX_MvmtBan_75_Epi_med <- IP_VAX_MvmtBan_75_Epi[,c(1,103)]
IP_VAX_3km_MvmtBan_90_Epi_med <- IP_VAX_3km_MvmtBan_90_Epi[,c(1,103)]
IP_VAX_3km_MvmtBan_75_Epi_med <- IP_VAX_3km_MvmtBan_75_Epi[,c(1,103)]
IP_VAX_10km_MvmtBan_90_Epi_med <- IP_VAX_10km_MvmtBan_90_Epi[,c(1,103)]
IP_VAX_10km_MvmtBan_75_Epi_med <- IP_VAX_10km_MvmtBan_75_Epi[,c(1,103)]
base_ShipmentsOff_Epi_med <- base_ShipmentsOff_Epi[,c(1,103)]


## Calculate the Upper 2.5% for each county
base_Epi$upper <- apply(base_Epi[grep("run_", colnames(base_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_90_Epi$upper <- apply(IP_MvmtBan_90_Epi[grep("run_", colnames(IP_MvmtBan_90_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_75_Epi$upper <- apply(IP_MvmtBan_75_Epi[grep("run_", colnames(IP_MvmtBan_75_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_90_Epi$upper <- apply(IP_DC_MvmtBan_90_Epi[grep("run_", colnames(IP_DC_MvmtBan_90_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_75_Epi$upper <- apply(IP_DC_MvmtBan_75_Epi[grep("run_", colnames(IP_DC_MvmtBan_75_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_90_Epi$upper <- apply(IP_VAX_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_MvmtBan_90_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_75_Epi$upper <- apply(IP_VAX_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_MvmtBan_75_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_Epi$upper <- apply(IP_VAX_3km_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_Epi$upper <- apply(IP_VAX_3km_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_Epi$upper <- apply(IP_VAX_10km_MvmtBan_90_Epi[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_Epi$upper <- apply(IP_VAX_10km_MvmtBan_75_Epi[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)
base_ShipmentsOff_Epi$upper <- apply(base_ShipmentsOff_Epi[grep("run_", colnames(base_ShipmentsOff_Epi))], 1, quantile, probs=0.975, na.rm = TRUE)

# Create a df with just the fips code and upper to map
base_Epi_up <- base_Epi[,c(1,104)]
IP_MvmtBan_90_Epi_up <- IP_MvmtBan_90_Epi[,c(1,104)]
IP_MvmtBan_75_Epi_up <- IP_MvmtBan_75_Epi[,c(1,104)]
IP_DC_MvmtBan_90_Epi_up <- IP_DC_MvmtBan_90_Epi[,c(1,104)]
IP_DC_MvmtBan_75_Epi_up <- IP_DC_MvmtBan_75_Epi[,c(1,104)]
IP_VAX_MvmtBan_90_Epi_up <- IP_VAX_MvmtBan_90_Epi[,c(1,104)]
IP_VAX_MvmtBan_75_Epi_up <- IP_VAX_MvmtBan_75_Epi[,c(1,104)]
IP_VAX_3km_MvmtBan_90_Epi_up <- IP_VAX_3km_MvmtBan_90_Epi[,c(1,104)]
IP_VAX_3km_MvmtBan_75_Epi_up <- IP_VAX_3km_MvmtBan_75_Epi[,c(1,104)]
IP_VAX_10km_MvmtBan_90_Epi_up <- IP_VAX_10km_MvmtBan_90_Epi[,c(1,104)]
IP_VAX_10km_MvmtBan_75_Epi_up <- IP_VAX_10km_MvmtBan_75_Epi[,c(1,104)]
base_ShipmentsOff_Epi_up <- base_ShipmentsOff_Epi[,c(1,104)]

# Create a vector for scale
Epi_scale <- c(base_Epi$upper, IP_MvmtBan_90_Epi$upper, IP_MvmtBan_75_Epi$upper, IP_DC_MvmtBan_90_Epi$upper, 
              IP_DC_MvmtBan_75_Epi$upper, IP_VAX_MvmtBan_90_Epi$upper, IP_VAX_MvmtBan_75_Epi$upper, 
              IP_VAX_3km_MvmtBan_90_Epi$upper, IP_VAX_3km_MvmtBan_75_Epi$upper, IP_VAX_10km_MvmtBan_90_Epi$upper, 
              IP_VAX_10km_MvmtBan_75_Epi$upper, base_ShipmentsOff_Epi$upper)
Epi_scale <- round(quantile(Epi_scale, probs = c(0, .21, .25, .3, .4, .5, .6, .7, 0.75, .8, 0.85, .9, 0.95, 0.99, 1), na.rm = TRUE), 0)
Epi_scale <- c(0, 1, 2, 4, 5, 10, 100, 500, 1000, 1200, 1400, 1535)

# Median maps
jpeg("Epi_Med_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Med_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_Epi_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

# Upper 2.5% maps
jpeg("Epi_Up_Base_min_Random.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Epi_Up_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_Epi_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Epi_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

#### Type of Spread (Detail Files) ####

detail.fnames <- list.files(path = path, recursive = TRUE, pattern = "_detail.txt", full.names = FALSE)
detail.fnames.90 <- detail.fnames[grep("90|75|NoCtrl|noctrl|base", detail.fnames)]

# Extract the infection route type out of the detail file
detail.fname.90 <- detail.fnames.90[1]
detail.res <- as.data.frame(fread(detail.fname.90, header = TRUE, select = c("SourceCounty", "InfRoute", "ControlPrevented")))
detail.res <- detail.res[detail.res$ControlPrevented == "none",]
detail.res <- as.data.frame(table(detail.res$SourceCounty, detail.res$InfRoute))
county.ls <- dcast(detail.res, Var1~Var2, value.var = "Freq")
colnames(county.ls) <- c("SourceCounty", "Local", "Ship")

county.ls$Type <- unlist(strsplit(detail.fname.90, "_flaps"))[1]
county.ls$Type <- unlist(strsplit(county.ls$Type, "/"))[2]

for(i in 2:length(detail.fnames.90)){
  detail.fname.90 <- detail.fnames.90[i]
  detail.res <- fread(detail.fname.90, header = TRUE, select = c("SourceCounty", "InfRoute", "ControlPrevented"))
  detail.res <- detail.res[detail.res$ControlPrevented == "none",]
  detail.res <- as.data.frame(table(detail.res$SourceCounty, detail.res$InfRoute))
  detail.res <- dcast(detail.res, Var1~Var2, value.var = "Freq")
  colnames(detail.res) <- c("SourceCounty", "Local", "Ship")
  detail.res$Type <- unlist(strsplit(detail.fname.90, "_flaps"))[1]
  detail.res$Type <- unlist(strsplit(detail.res$Type, "/"))[2]
  detail.res <- as.data.frame(detail.res)
  county.ls <- merge(county.ls, detail.res, by = "SourceCounty", all=TRUE)
  county.ls[is.na(county.ls)] <- 0
  
}

names(county.ls)[seq(4,ncol(county.ls),3)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.ls))/2))))
names(county.ls)[seq(2,ncol(county.ls),3)] <- paste0("local_", sprintf("%03d", seq(1:((ncol(county.ls))/2))))
names(county.ls)[seq(3,ncol(county.ls),3)] <- paste0("ship_", sprintf("%03d", seq(1:((ncol(county.ls))/2))))

write.csv(county.ls, "county_ls.csv")
# county.ls <- read.csv("county_ls.csv", header = TRUE)
# county.ls <- county.ls[,-1]
base_ls <- new_base_ls

## Create individual data frames for each run type and local spread
df <- county.ls[,c(1,2:301)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,302:601)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,602:901)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,902:1201)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,1202:1501)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,1502:1801)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,1802:2101)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,2102:2401)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,2402:2701)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,2702:3001)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,3002:3301)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,3302:3601)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

df <- county.ls[,c(1,3602:3901)]
df1 <- df[df$type_1201 != 0 & !is.na(df$type_1201),]
newname <- df1$type_1201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_ls"), df.new)

## calculate sum of each run's local spread
base_ls$local <- apply(base_ls[grep("local_", colnames(base_ls))], 1, sum, na.rm = TRUE)
IP_MvmtBan_90_ls$local <- apply(IP_MvmtBan_90_ls[grep("local_", colnames(IP_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_MvmtBan_75_ls$local <- apply(IP_MvmtBan_75_ls[grep("local_", colnames(IP_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_DC_MvmtBan_90_ls$local <- apply(IP_DC_MvmtBan_90_ls[grep("local_", colnames(IP_DC_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_DC_MvmtBan_75_ls$local <- apply(IP_DC_MvmtBan_75_ls[grep("local_", colnames(IP_DC_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_MvmtBan_90_ls$local <- apply(IP_VAX_MvmtBan_90_ls[grep("local_", colnames(IP_VAX_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_MvmtBan_75_ls$local <- apply(IP_VAX_MvmtBan_75_ls[grep("local_", colnames(IP_VAX_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_ls$local <- apply(IP_VAX_3km_MvmtBan_90_ls[grep("local_", colnames(IP_VAX_3km_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_ls$local <- apply(IP_VAX_3km_MvmtBan_75_ls[grep("local_", colnames(IP_VAX_3km_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_ls$local <- apply(IP_VAX_10km_MvmtBan_90_ls[grep("local_", colnames(IP_VAX_10km_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_ls$local <- apply(IP_VAX_10km_MvmtBan_75_ls[grep("local_", colnames(IP_VAX_10km_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
base_ShipmentsOff_ls$local <- apply(base_ShipmentsOff_ls[grep("local_", colnames(base_ShipmentsOff_ls))], 1, sum, na.rm = TRUE)

## calculate sum of each run's shipment spread
base_ls$ship <- apply(base_ls[grep("ship_", colnames(base_ls))], 1, sum, na.rm = TRUE)
IP_MvmtBan_90_ls$ship <- apply(IP_MvmtBan_90_ls[grep("ship_", colnames(IP_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_MvmtBan_75_ls$ship <- apply(IP_MvmtBan_75_ls[grep("ship_", colnames(IP_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_DC_MvmtBan_90_ls$ship <- apply(IP_DC_MvmtBan_90_ls[grep("ship_", colnames(IP_DC_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_DC_MvmtBan_75_ls$ship <- apply(IP_DC_MvmtBan_75_ls[grep("ship_", colnames(IP_DC_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_MvmtBan_90_ls$ship <- apply(IP_VAX_MvmtBan_90_ls[grep("ship_", colnames(IP_VAX_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_MvmtBan_75_ls$ship <- apply(IP_VAX_MvmtBan_75_ls[grep("ship_", colnames(IP_VAX_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_ls$ship <- apply(IP_VAX_3km_MvmtBan_90_ls[grep("ship_", colnames(IP_VAX_3km_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_ls$ship <- apply(IP_VAX_3km_MvmtBan_75_ls[grep("ship_", colnames(IP_VAX_3km_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_ls$ship <- apply(IP_VAX_10km_MvmtBan_90_ls[grep("ship_", colnames(IP_VAX_10km_MvmtBan_90_ls))], 1, sum, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_ls$ship <- apply(IP_VAX_10km_MvmtBan_75_ls[grep("ship_", colnames(IP_VAX_10km_MvmtBan_75_ls))], 1, sum, na.rm = TRUE)
base_ShipmentsOff_ls$ship <- apply(base_ShipmentsOff_ls[grep("ship_", colnames(base_ShipmentsOff_ls))], 1, sum, na.rm = TRUE)

# Add in the proportion of spread that is local
base_ls$localProp <- round((base_ls$local/(base_ls$local + base_ls$ship))*100,0)
IP_MvmtBan_90_ls$localProp <- round((IP_MvmtBan_90_ls$local/(IP_MvmtBan_90_ls$local + IP_MvmtBan_90_ls$ship))*100,0)
IP_MvmtBan_75_ls$localProp <- round((IP_MvmtBan_75_ls$local/(IP_MvmtBan_75_ls$local + IP_MvmtBan_75_ls$ship))*100,0)
IP_DC_MvmtBan_90_ls$localProp <- round((IP_DC_MvmtBan_90_ls$local/(IP_DC_MvmtBan_90_ls$local + IP_DC_MvmtBan_90_ls$ship))*100,0)
IP_DC_MvmtBan_75_ls$localProp <- round((IP_DC_MvmtBan_75_ls$local/(IP_DC_MvmtBan_75_ls$local + IP_DC_MvmtBan_75_ls$ship))*100,0)
IP_VAX_MvmtBan_90_ls$localProp <- round((IP_VAX_MvmtBan_90_ls$local/(IP_VAX_MvmtBan_90_ls$local + IP_VAX_MvmtBan_90_ls$ship))*100,0)
IP_VAX_MvmtBan_75_ls$localProp <- round((IP_VAX_MvmtBan_90_ls$local/(IP_VAX_MvmtBan_90_ls$local + IP_VAX_MvmtBan_90_ls$ship))*100,0)
IP_VAX_3km_MvmtBan_90_ls$localProp <- round((IP_VAX_3km_MvmtBan_90_ls$local/(IP_VAX_3km_MvmtBan_90_ls$local + IP_VAX_3km_MvmtBan_90_ls$ship))*100,0)
IP_VAX_3km_MvmtBan_75_ls$localProp <- round((IP_VAX_3km_MvmtBan_90_ls$local/(IP_VAX_3km_MvmtBan_90_ls$local + IP_VAX_3km_MvmtBan_90_ls$ship))*100,0)
IP_VAX_10km_MvmtBan_90_ls$localProp <- round((IP_VAX_10km_MvmtBan_90_ls$local/(IP_VAX_10km_MvmtBan_90_ls$local + IP_VAX_10km_MvmtBan_90_ls$ship))*100,0)
IP_VAX_10km_MvmtBan_75_ls$localProp <- round((IP_VAX_10km_MvmtBan_90_ls$local/(IP_VAX_10km_MvmtBan_90_ls$local + IP_VAX_10km_MvmtBan_90_ls$ship))*100,0)
base_ShipmentsOff_ls$localProp <- round((base_ShipmentsOff_ls$local/(base_ShipmentsOff_ls$local + base_ShipmentsOff_ls$ship))*100,0)

# Create a data.frame of just two columns for mapping
base_local <- base_ls[,c("SourceCounty", "localProp")]
IP_MvmtBan_90_local <- IP_MvmtBan_90_ls[,c("SourceCounty", "localProp")]
IP_MvmtBan_75_local <- IP_MvmtBan_75_ls[,c("SourceCounty", "localProp")]
IP_DC_MvmtBan_90_local <- IP_DC_MvmtBan_90_ls[,c("SourceCounty", "localProp")]
IP_DC_MvmtBan_75_local <- IP_DC_MvmtBan_75_ls[,c("SourceCounty", "localProp")]
IP_VAX_MvmtBan_90_local <- IP_VAX_MvmtBan_90_ls[,c("SourceCounty", "localProp")]
IP_VAX_MvmtBan_75_local <- IP_VAX_MvmtBan_75_ls[,c("SourceCounty", "localProp")]
IP_VAX_3km_MvmtBan_90_local <- IP_VAX_3km_MvmtBan_90_ls[,c("SourceCounty", "localProp")]
IP_VAX_3km_MvmtBan_75_local <- IP_VAX_3km_MvmtBan_75_ls[,c("SourceCounty", "localProp")]
IP_VAX_10km_MvmtBan_90_local <- IP_VAX_10km_MvmtBan_90_ls[,c("SourceCounty", "localProp")]
IP_VAX_10km_MvmtBan_75_local <- IP_VAX_10km_MvmtBan_75_ls[,c("SourceCounty", "localProp")]
base_ShipmentsOff_local <- base_ShipmentsOff_ls[,c("SourceCounty", "localProp")]

# Create a vector for scale
local_scale <- round(c(0, 25, 50, 60, 70, 80, 85, 90, 95, 97.5, 99, 100),2)

# Create maps of Local spreads
jpeg("Local_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Local_ShipOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_local, county.border.col = NA, state.border.col = "gray20", missing.include = TRUE, color.break.type = "values", 
            color.break.values = local_scale, color.sequence = color_bluepurple, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()


# Create maps of Shipment spreads
jpeg("Ship_Base.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Ship_IP90.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Ship_IPDC90.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Ship_ShipOff.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Ship_IP75.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Ship_IPDC75.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_ship, county.border.col = NA, state.border.col = "gray40", missing.include = TRUE, color.break.type = "values", 
            color.break.values = scale_100, color.sequence = color_orange, legend.spacing = 5.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

#### Risk (Detail Files) ####

# County Risk

# read the summary file
sum.file <- as.data.frame(fread(summary.fnames.90[1], header = TRUE, select = c("Rep", "Seed_FIPS")))

# read and process the detail file 
det.file <- as.data.frame(fread(detail.fnames.90[1], header = TRUE, select = c("Rep","SourceCounty","ExposedCounty", "ControlPrevented")))
realized <- det.file[det.file$ControlPrevented == "none",] # successful exposures
realized2 <- unique(as.data.table(realized),by=c("Rep","ExposedCounty")) # only count a county once per rep
realized3 <- realized2[which(realized2$SourceCounty!=realized2$ExposedCounty),] # do not include the seed county
# Note though, if there's a situation where county 1 (the seed) infects co. 2, which infects co. 3, which 
# infects co. 1, 1 would be counted. The below two lines fix that case. 
realized4=merge(sum.file,realized3,by="Rep",all=T)
realized5=realized4[realized4$Seed_FIPS!=realized4$ExposedCounty,]

fr= as.data.frame(table(realized5$ExposedCounty))
fr$Type <- unlist(strsplit(summary.fnames.90[1], "_flaps"))[1]
fr$Type <- unlist(strsplit(fr$Type, "/"))[2]

for (file in 2:length(detail.fnames.90)){
  
  sum.file <- as.data.frame(fread(summary.fnames.90[file], header = TRUE, select = c("Rep", "Seed_FIPS")))
  
  # read and process the detail file 
  det.file <- as.data.frame(fread(detail.fnames.90[file], header = TRUE, select = c("Rep","SourceCounty","ExposedCounty", "ControlPrevented")))
  realized <- det.file[det.file$ControlPrevented == "none",] 
  realized2 <- unique(as.data.table(realized),by=c("Rep","ExposedCounty")) 
  realized3 <- realized2[which(realized2$SourceCounty!=realized2$ExposedCounty),] 
  realized4=merge(sum.file,realized3,by="Rep",all=T)
  realized5=realized4[realized4$Seed_FIPS!=realized4$ExposedCounty,]
  f= as.data.frame(table(realized5$ExposedCounty))
  
  f$Type <- unlist(strsplit(summary.fnames.90[file], "_flaps"))[1]
  f$Type <- unlist(strsplit(f$Type, "/"))[2]
  
  fr=merge(fr,f, by="Var1",all=T)
  print(file)
}

# Formatting
names(fr)[1]="fips"
fr=merge(county.fips,fr,by="fips",all=T)
fr <- unique(as.data.table(fr),by=c("fips"))

fr=as.data.frame(fr)

names(fr)[seq(3,ncol(fr),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(fr))/2-1))))
names(fr)[seq(4,ncol(fr),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(fr))/2-1))))

fr$fips[which(fr$fips==46113)]<- 46102 # A county name/code change 

fr[, seq(3,ncol(fr),2)][is.na(fr[, seq(3,ncol(fr),2)])] <- 0
# replacing missing run types
fr=na.locf(fr,na.rm=F) 
fr=na.locf(fr,na.rm=F,fromLast = T)

# To make this into a proportion divide by 3049 (number of counties)
fr[,seq(3,ncol(fr),2)] <- fr[,seq(3,ncol(fr),2)]/3049

###
write.csv(fr, "CountyRisk_Binary.csv",row.names=F)
fr <- read.csv("CountyRisk_Binary.csv", header = TRUE)

## Create individual data frames for each run type

df <- fr[,c(1,2,3:202)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,203:402)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,403:602)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,603:802)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,803:1002)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,1003:1202)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,1203:1402)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,1403:1602)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,1603:1802)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,1803:2002)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,2003:2202)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

df <- fr[,c(1,2,2203:2402)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_fr"), df.new)

## calculate median of each run
base_fr$median <- apply(base_fr[grep("run_", colnames(base_fr))], 1, median, na.rm = TRUE)
IP_MvmtBan_90_fr$median <- apply(IP_MvmtBan_90_fr[grep("run_", colnames(IP_MvmtBan_90_fr))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_90_fr$median <- apply(IP_DC_MvmtBan_90_fr[grep("run_", colnames(IP_DC_MvmtBan_90_fr))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_90_fr$median <- apply(IP_VAX_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_MvmtBan_90_fr))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_fr$median <- apply(IP_VAX_3km_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_fr))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_fr$median <- apply(IP_VAX_10km_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_fr))], 1, median, na.rm = TRUE)
IP_MvmtBan_75_fr$median <- apply(IP_MvmtBan_75_fr[grep("run_", colnames(IP_MvmtBan_75_fr))], 1, median, na.rm = TRUE)
IP_DC_MvmtBan_75_fr$median <- apply(IP_DC_MvmtBan_75_fr[grep("run_", colnames(IP_DC_MvmtBan_75_fr))], 1, median, na.rm = TRUE)
IP_VAX_MvmtBan_75_fr$median <- apply(IP_VAX_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_MvmtBan_75_fr))], 1, median, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_fr$median <- apply(IP_VAX_3km_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_fr))], 1, median, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_fr$median <- apply(IP_VAX_10km_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_fr))], 1, median, na.rm = TRUE)
base_ShipmentsOff_fr$median <- apply(base_ShipmentsOff_fr[grep("run_", colnames(base_ShipmentsOff_fr))], 1, median, na.rm = TRUE)


# Create a data.frame with just the fips code and median to map
base_fr_med <- base_fr[,c(1,103)]
IP_MvmtBan_90_fr_med <- IP_MvmtBan_90_fr[,c(1,103)]
IP_DC_MvmtBan_90_fr_med <- IP_DC_MvmtBan_90_fr[,c(1,103)]
IP_VAX_MvmtBan_90_fr_med <- IP_VAX_MvmtBan_90_fr[,c(1,103)]
IP_VAX_3km_MvmtBan_90_fr_med <- IP_VAX_3km_MvmtBan_90_fr[,c(1,103)]
IP_VAX_10km_MvmtBan_90_fr_med <- IP_VAX_10km_MvmtBan_90_fr[,c(1,103)]
IP_MvmtBan_75_fr_med <- IP_MvmtBan_75_fr[,c(1,103)]
IP_DC_MvmtBan_75_fr_med <- IP_DC_MvmtBan_75_fr[,c(1,103)]
IP_VAX_MvmtBan_75_fr_med <- IP_VAX_MvmtBan_75_fr[,c(1,103)]
IP_VAX_3km_MvmtBan_75_fr_med <- IP_VAX_3km_MvmtBan_75_fr[,c(1,103)]
IP_VAX_10km_MvmtBan_75_fr_med <- IP_VAX_10km_MvmtBan_75_fr[,c(1,103)]
base_ShipmentsOff_fr_med <- base_ShipmentsOff_fr[,c(1,103)]

## Calculate the Upper 2.5% for each county
base_fr$upper <- apply(base_fr[grep("run_", colnames(base_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_90_fr$upper <- apply(IP_MvmtBan_90_fr[grep("run_", colnames(IP_MvmtBan_90_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_90_fr$upper <- apply(IP_DC_MvmtBan_90_fr[grep("run_", colnames(IP_DC_MvmtBan_90_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_90_fr$upper <- apply(IP_VAX_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_MvmtBan_90_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_90_fr$upper <- apply(IP_VAX_3km_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_90_fr$upper <- apply(IP_VAX_10km_MvmtBan_90_fr[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_MvmtBan_75_fr$upper <- apply(IP_MvmtBan_75_fr[grep("run_", colnames(IP_MvmtBan_75_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_DC_MvmtBan_75_fr$upper <- apply(IP_DC_MvmtBan_75_fr[grep("run_", colnames(IP_DC_MvmtBan_75_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_MvmtBan_75_fr$upper <- apply(IP_VAX_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_MvmtBan_75_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_3km_MvmtBan_75_fr$upper <- apply(IP_VAX_3km_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
IP_VAX_10km_MvmtBan_75_fr$upper <- apply(IP_VAX_10km_MvmtBan_75_fr[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_fr))], 1, quantile, probs=0.975, na.rm = TRUE)
base_ShipmentsOff_fr$upper <- apply(base_ShipmentsOff_fr[grep("run_", colnames(base_ShipmentsOff_fr))], 1, quantile, probs=0.975, na.rm = TRUE)

# Create a data.frame with just the fips code and upper to map
base_fr_up <- base_fr[,c(1,104)]
IP_MvmtBan_90_fr_up <- IP_MvmtBan_90_fr[,c(1,104)]
IP_DC_MvmtBan_90_fr_up <- IP_DC_MvmtBan_90_fr[,c(1,104)]
IP_VAX_MvmtBan_90_fr_up <- IP_VAX_MvmtBan_90_fr[,c(1,104)]
IP_VAX_3km_MvmtBan_90_fr_up <- IP_VAX_3km_MvmtBan_90_fr[,c(1,104)]
IP_VAX_10km_MvmtBan_90_fr_up <- IP_VAX_10km_MvmtBan_90_fr[,c(1,104)]
IP_MvmtBan_75_fr_up <- IP_MvmtBan_75_fr[,c(1,104)]
IP_DC_MvmtBan_75_fr_up <- IP_DC_MvmtBan_75_fr[,c(1,104)]
IP_VAX_MvmtBan_75_fr_up <- IP_VAX_MvmtBan_75_fr[,c(1,104)]
IP_VAX_3km_MvmtBan_75_fr_up <- IP_VAX_3km_MvmtBan_75_fr[,c(1,104)]
IP_VAX_10km_MvmtBan_75_fr_up <- IP_VAX_10km_MvmtBan_75_fr[,c(1,104)]
base_ShipmentsOff_fr_up <- base_ShipmentsOff_fr[,c(1,104)]

# Create a vector for scale
Fr_uppers <- c(base_fr$upper, IP_MvmtBan_90_fr$upper, IP_DC_MvmtBan_90_fr$upper, IP_VAX_MvmtBan_90_fr$upper, 
               IP_VAX_3km_MvmtBan_90_fr$upper, IP_VAX_10km_MvmtBan_90_fr$upper)
Fr_upper_values <- round(quantile(Fr_uppers, probs = c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1), na.rm = TRUE), 3)
Fr_upper_values <- c(0.000, 0.002, 0.007, 0.011, 0.020, 0.033, 0.036, 0.041)
Fr_upper_values <- sort(unique(Fr_upper_values))
Fr_upper_values[length(Fr_upper_values)]=Fr_upper_values[length(Fr_upper_values)]+0.001 # to avoid rounding to 3 decimal places issues with 


# Median maps
jpeg("Fr_Med_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5,
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Med_Base_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_fr_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 6.5,
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()


# Upper 2.5% maps
jpeg("Fr_Up_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_fr_up,county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_fr_up,county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_fr_up,county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_fr_up,county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()

jpeg("Fr_Up_Base_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_fr_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Fr_upper_values, color.sequence = color_YlOrRd, legend.spacing = 4.5, 
            legend.shrink = 0.5, legend.width = 1,legend.digits=3)
dev.off()


#### DC Culls implemented (Summary Files) ####

county.dc <- county.fips

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  if(!grepl("DC",summary.fname.90)){next} else {
    summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "cullImplementedDCSubset"))
    summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
    summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
    summary.res <- as.data.frame(summary.res)
    county.dc <- merge(county.dc, summary.res, by.x = "fips", by.y = "Seed_FIPS", all=TRUE)
  }
}
dup_fips <- rownames(county.dc)[duplicated(county.dc$fips) == TRUE]
dup_name <- rownames(county.dc)[duplicated(county.dc$polyname) == TRUE]
county.dc <- county.dc[!(rownames(county.dc) %in% dup_fips),]
county.dc <- county.dc[!(rownames(county.dc) %in% dup_name),]

names(county.dc)[seq(4,ncol(county.dc),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.dc))/2-1))))
names(county.dc)[seq(3,ncol(county.dc),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.dc))/2-1))))

write.csv(county.dc, "county_dc.csv")
# county.dc <- read.csv("county_dc.csv", header = TRUE)
# county.dc <- county.dc[,-1]

df <- county.dc[,c(1,2,3:202)]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(df[1,3], "_DCcull"), df.new)

# Convert to long format for plotting
long.county.dc <- reshape(county.dc, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                          varying = list(c(grep("type_", colnames(county.dc))), c(grep("run_", colnames(county.dc)))))
rownames(long.county.dc) <- seq(1:nrow(long.county.dc))
long.county.dc <- long.county.dc[,-3]
long.county.dc <- long.county.dc[!is.na(long.county.dc$Value),]

# Rename the types of runs
long.county.dc$type <- as.factor(long.county.dc$type)
levels(long.county.dc$type)[levels(long.county.dc$type) == "IP_DC_MvmtBan_90_"] <- "IP & DC Cull"

head(long.county.dc)

jpeg("DC_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dc, aes(x = long.county.dc$type, y = long.county.dc$Value, color = long.county.dc$type)) + 
  scale_x_discrete(limits = c("IP & DC Cull")) + labs(x = NULL, y = "DC Culls", title = expression("DC Culls")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

# Divide the data
table(long.county.dc$Value)
jpeg("DC_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.dc$Value[long.county.dc$Value > 2000], breaks = 25, xlab = "DC Culls", main = expression("DC Culls" > "2000 Subset"))
dev.off()
jpeg("DC_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.dc$Value[long.county.dc$Value <= 2000], breaks = 25, xlab = "DC Culls", main = expression("DC Culls" <= "2000 Subset"))
dev.off()

long.county.dc.less5 <- long.county.dc[long.county.dc$Value > 5,]
long.county.dc.low <- long.county.dc[long.county.dc$Value <= 2000,]
long.county.dc.high <- long.county.dc[long.county.dc$Value > 2000,]

jpeg("DC_Less5_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dc.less5, aes(x = long.county.dc.less5$type, y = long.county.dc.less5$Value, color = long.county.dc.less5$type)) + 
  scale_x_discrete(limits = c("IP & DC Cull")) + labs(x = NULL, y = "DC Culls", title = expression("DC Culls" > "5 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("DC_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dc.low, aes(x = long.county.dc.low$type, y = long.county.dc.low$Value, color = long.county.dc.low$type)) + 
  scale_x_discrete(limits = c("IP & DC Cull")) + labs(x = NULL, y = "DC Culls", title = expression("DC Culls" <= "2000 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("DC_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.dc.high, aes(x = long.county.dc.high$type, y = long.county.dc.high$Value, color = long.county.dc.high$type)) + 
  scale_x_discrete(limits = c("IP & DC Cull")) + labs(x = NULL, y = "DC Culls", title = expression("DC Culls" > "2000 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

# Get the medians for each control type

IP_DC90_DC_Low <- median(long.county.dc.low$Value[long.county.dc.low$type == "IP & DC Cull"])
IP_DC90_DC_High <- median(long.county.dc.high$Value[long.county.dc.high$type == "IP & DC Cull"])

DC_Medians_90_Low <- rbind("NA", "NA", IP_DC90_DC_Low, "NA", "NA", "NA")
DC_Medians_90_High <- rbind("NA", "NA", IP_DC90_DC_High, "NA", "NA", "NA")

#### IP Culls implemented (Summary Files) ####

county.ip <- county.fips

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "cullImplemented"))
  summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
  summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
  summary.res <- as.data.frame(summary.res)
  summary.res$cullImplemented[summary.res$cullImplemented == 0] <- NA
  county.ip <- merge(county.ip, summary.res, by.x = "fips", by.y = "Seed_FIPS", all=TRUE)
}
dup_fips <- rownames(county.ip)[duplicated(county.ip$fips) == TRUE]
dup_name <- rownames(county.ip)[duplicated(county.ip$polyname) == TRUE]
county.ip <- county.ip[!(rownames(county.ip) %in% dup_fips),]
county.ip <- county.ip[!(rownames(county.ip) %in% dup_name),]

names(county.ip)[seq(4,ncol(county.ip),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.ip))/2-1))))
names(county.ip)[seq(3,ncol(county.ip),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.ip))/2-1))))

write.csv(county.ip, "county_ip.csv")
# county.ip <- read.csv("county_ip.csv", header = TRUE)
# county.ip <- county.ip[,-1]

# Convert to long format for plotting
long.county.ip <- reshape(county.ip, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                          varying = list(c(grep("type_", colnames(county.ip))), c(grep("run_", colnames(county.ip)))))
rownames(long.county.ip) <- seq(1:nrow(long.county.ip))
long.county.ip <- long.county.ip[,-3]
long.county.ip <- long.county.ip[!is.na(long.county.ip$Value),]

# Rename the types of runs
long.county.ip$type <- as.factor(long.county.ip$type)
levels(long.county.ip$type)[levels(long.county.ip$type) == "IP_MvmtBan_90_"] <- "IP Cull"
levels(long.county.ip$type)[levels(long.county.ip$type) == "IP_DC_MvmtBan_90_"] <- "IP & DC Cull"
levels(long.county.ip$type)[levels(long.county.ip$type) == "base_ShipmentsOff"] <- "IP Cull & DC Vax"
levels(long.county.ip$type)[levels(long.county.ip$type) == "IP_MvmtBan_75"] <- "IP Cull & 3km Vax"
levels(long.county.ip$type)[levels(long.county.ip$type) == "IP_DC_MvmtBan_75"] <- "IP Cull & 10km Vax"

head(long.county.ip)

# Divide the data
long.county.ip.less10 <- long.county.ip[long.county.ip$Value > 10,]
long.county.ip.low <- long.county.ip[long.county.ip$Value <= 3000,]
long.county.ip.high <- long.county.ip[long.county.ip$Value > 3000,]

table(long.county.ip$Value)
jpeg("IP_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.ip$Value[long.county.ip$Value > 3000], breaks = 25, xlab = "DC Culls", main = expression("IP Culls" > "3000 Subset"))
dev.off()
jpeg("IP_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.ip$Value[long.county.ip$Value <= 3000], breaks = 25, xlab = "DC Culls", main = expression("IP Culls" <= "3000 Subset"))
dev.off()

jpeg("IP_Less10_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.ip.less10, aes(x = long.county.ip.less10$type, y = long.county.ip.less10$Value, color = long.county.ip.less10$type)) + 
  scale_x_discrete(limits = c("IP Cull", "IP & DC Cull", "IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "IP Culls", title = expression("IPs Culled" > "10 Subset")) + scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("IP_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.ip.low, aes(x = long.county.ip.low$type, y = long.county.ip.low$Value, color = long.county.ip.low$type)) + 
  scale_x_discrete(limits = c("IP Cull", "IP & DC Cull", "IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "IP Culls", title = expression("IPs Culled" <= "3000 Subset")) + scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("IP_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.ip.high, aes(x = long.county.ip.high$type, y = long.county.ip.high$Value, color = long.county.ip.high$type)) + 
  scale_x_discrete(limits = c("IP Cull", "IP & DC Cull", "IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "IP Culls", title = expression("IPs Culled" > "3000 Subset")) + scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(legend.position = "none")
dev.off()

# Calculate the medians for each control type
IP90_IP_Low <- median(long.county.ip.low$Value[long.county.ip.low$type == "IP Cull"])
IP90_IP_High <- median(long.county.ip.high$Value[long.county.ip.high$type == "IP Cull"])
IP_DC90_IP_Low <- median(long.county.ip.low$Value[long.county.ip.low$type == "IP & DC Cull"]) - IP_DC90_DC_Low
IP_DC90_IP_High <- median(long.county.ip.high$Value[long.county.ip.high$type == "IP & DC Cull"]) - IP_DC90_DC_High
IPVAX90_IP_Low <- median(long.county.ip.low$Value[long.county.ip.low$type == "IP Cull & DC Vax"])
IPVAX90_IP_High <- median(long.county.ip.high$Value[long.county.ip.high$type == "IP Cull & DC Vax"])
IPVAX90_3_IP_Low <- median(long.county.ip.low$Value[long.county.ip.low$type == "IP Cull & 3km Vax"])
IPVAX90_3_IP_High <- median(long.county.ip.high$Value[long.county.ip.high$type == "IP Cull & 3km Vax"])
IPVAX90_10_IP_Low <- median(long.county.ip.low$Value[long.county.ip.low$type == "IP Cull & 10km Vax"])
IPVAX90_10_IP_High <- median(long.county.ip.high$Value[long.county.ip.high$type == "IP Cull & 10km Vax"])

IP_Medians_90_Low <- rbind("NA", IP90_IP_Low, IP_DC90_IP_Low, IPVAX90_IP_Low, IPVAX90_3_IP_Low, IPVAX90_10_IP_Low)
IP_Medians_90_High <- rbind("NA", IP90_IP_High, IP_DC90_IP_High, IPVAX90_IP_High, IPVAX90_3_IP_High, IPVAX90_10_IP_High)


#### Vaccines implemented (Summary Files) ####

county.vax <- county.fips

for(i in 1:length(summary.fnames.90)){
  summary.fname.90 <- summary.fnames.90[i]
  if(!grepl("VAX",summary.fname.90))next
  summary.res <- fread(summary.fname.90, header = TRUE, select = c("Seed_FIPS", "vaxEffective"))
  summary.res$Type <- unlist(strsplit(summary.fname.90, "_flaps"))[1]
  summary.res$Type <- unlist(strsplit(summary.res$Type, "/"))[2]
  summary.res <- as.data.frame(summary.res)
  county.vax <- merge(county.vax, summary.res, by.x = "fips", by.y = "Seed_FIPS", all=TRUE)
}
dup_fips <- rownames(county.vax)[duplicated(county.vax$fips) == TRUE]
dup_name <- rownames(county.vax)[duplicated(county.vax$polyname) == TRUE]
county.vax <- county.vax[!(rownames(county.vax) %in% dup_fips),]
county.vax <- county.vax[!(rownames(county.vax) %in% dup_name),]

names(county.vax)[seq(4,ncol(county.vax),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.vax))/2-1))))
names(county.vax)[seq(3,ncol(county.vax),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.vax))/2-1))))

write.csv(county.vax, "county_vax.csv")
county.vax <- read.csv("county_vax.csv", header = TRUE)
county.vax <- county.vax[,-1]

# Convert to long format for plotting
length(unique(county.vax$fips))
length(county.vax$fips)
dup_fips <- rownames(county.vax)[duplicated(county.vax$fips) == TRUE]
dup_name <- rownames(county.vax)[duplicated(county.vax$polyname) == TRUE]
county.vax <- county.vax[!(rownames(county.vax) %in% dup_fips),]
county.vax <- county.vax[!(rownames(county.vax) %in% dup_name),]

long.county.vax <- reshape(county.vax, idvar = c("fips", "polyname"), direction = "long", v.names = c("type", "Value"), 
                           varying = list(c(grep("type_", colnames(county.vax))), c(grep("run_", colnames(county.vax)))))
rownames(long.county.vax) <- seq(1:nrow(long.county.vax))
long.county.vax <- long.county.vax[,-3]
long.county.vax <- long.county.vax[!is.na(long.county.vax$Value),]

# Rename the types of runs
long.county.vax$type <- as.factor(long.county.vax$type)
levels(long.county.vax$type)[levels(long.county.vax$type) == "base_ShipmentsOff"] <- "IP Cull & DC Vax"
levels(long.county.vax$type)[levels(long.county.vax$type) == "IP_MvmtBan_75"] <- "IP Cull & 3km Vax"
levels(long.county.vax$type)[levels(long.county.vax$type) == "IP_DC_MvmtBan_75"] <- "IP Cull & 10km Vax"

head(long.county.vax)

jpeg("Vax_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.vax, aes(x = long.county.vax$type, y = long.county.vax$Value, color = long.county.vax$type)) + 
  scale_x_discrete(limits = c("IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + labs(x = NULL, y = "Vaccinations", title = expression("Vaccinations")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

# Divide the data

table(long.county.vax$Value)
jpeg("Vax_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.vax$Value[long.county.vax$Value > 3000], breaks = 25, xlab = "Vaccinations", main = expression("Vaccinations" > "3000 Subset"))
dev.off()
jpeg("Vax_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
hist(long.county.vax$Value[long.county.vax$Value <= 3000], breaks = 25, xlab = "Vaccinations", main = expression("Vaccinations" <= "3000 Subset"))
dev.off()

long.county.vax.less10 <- long.county.vax[long.county.vax$Value > 10,]
long.county.vax.low <- long.county.vax[long.county.vax$Value <= 3000,]
long.county.vax.high <- long.county.vax[long.county.vax$Value > 3000,]

jpeg("Vax_Less10_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.vax.less10, aes(x = long.county.vax.less10$type, y = long.county.vax.less10$Value, color = long.county.vax.less10$type)) + 
  scale_x_discrete(limits = c("IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "Vaccinations", title = expression("Vaccinations" > "10 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("Vax_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.vax.low, aes(x = long.county.vax.low$type, y = long.county.vax.low$Value, color = long.county.vax.low$type)) + 
  scale_x_discrete(limits = c("IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "Vaccinations", title = expression("Vaccinations" <= "3000 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

jpeg("Vax_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.vax.high, aes(x = long.county.vax.high$type, y = long.county.vax.high$Value, color = long.county.vax.high$type)) + 
  scale_x_discrete(limits = c("IP Cull & DC Vax", "IP Cull & 3km Vax", "IP Cull & 10km Vax")) + 
  labs(x = NULL, y = "Vaccinations", title = expression("Vaccinations" > "3000 Subset")) +
  scale_color_manual(values = cbPalette) + theme_bw() + theme(legend.position = "none")
dev.off()

# Get the medians for each control type
IPVAX90_Vacc_Low <- median(long.county.vax.low$Value[long.county.vax.low$type == "IP Cull & DC Vax"])
IPVAX90_Vacc_High <- median(long.county.vax.high$Value[long.county.vax.high$type == "IP Cull & DC Vax"])
IPVAX90_3_Vacc_Low <- median(long.county.vax.low$Value[long.county.vax.low$type == "IP Cull & 3km Vax"])
IPVAX90_3_Vacc_High <- median(long.county.vax.high$Value[long.county.vax.high$type == "IP Cull & 3km Vax"])
IPVAX90_10_Vacc_Low <- median(long.county.vax.low$Value[long.county.vax.low$type == "IP Cull & 10km Vax"])
IPVAX90_10_Vacc_High <- median(long.county.vax.high$Value[long.county.vax.high$type == "IP Cull & 10km Vax"])

Vacc_Medians_90_Low <- rbind("NA", "NA", "NA", IPVAX90_Vacc_Low, IPVAX90_3_Vacc_Low, IPVAX90_10_Vacc_Low)
Vacc_Medians_90_High <- rbind("NA", "NA", "NA", IPVAX90_Vacc_High, IPVAX90_3_Vacc_High, IPVAX90_10_Vacc_High)

#### Number of Animals Infected (Detail Files) ####

# Read in and format all FLAPS files 
# min
f1_min <- fread("flaps12_min_0001.txt") 
f1_min$anim <- f1_min$V5 + f1_min$V6
f1_min$V2[which(f1_min$V2==46113)]<- 46102 # Fixing a county name/code change (Shannon/Lakota County SD)

f2_min <- fread("flaps12_min_0002.txt")
f2_min$anim <- f2_min$V5 + f2_min$V6
f2_min$V2[which(f2_min$V2==46113)]<- 46102

f3_min <- fread("flaps12_min_0003.txt")
f3_min$anim <- f3_min$V5 + f3_min$V6
f3_min$V2[which(f3_min$V2==46113)]<- 46102

f4_min <- fread("flaps12_min_0004.txt")
f4_min$anim <- f4_min$V5 + f4_min$V6
f4_min$V2[which(f4_min$V2==46113)]<- 46102

f5_min <- fread("flaps12_min_0005.txt")
f5_min$anim <- f5_min$V5 + f5_min$V6
f5_min$V2[which(f5_min$V2==46113)]<- 46102

f6_min <- fread("flaps12_min_0006.txt")
f6_min$anim <- f6_min$V5 + f6_min$V6
f6_min$V2[which(f6_min$V2==46113)]<- 46102

f7_min <- fread("flaps12_min_0007.txt")
f7_min$anim <- f7_min$V5 + f7_min$V6
f7_min$V2[which(f7_min$V2==46113)]<- 46102

f8_min <- fread("flaps12_min_0008.txt")
f8_min$anim <- f8_min$V5 + f8_min$V6
f8_min$V2[which(f8_min$V2==46113)]<- 46102

f9_min <- fread("flaps12_min_0009.txt")
f9_min$anim <- f9_min$V5 + f9_min$V6
f9_min$V2[which(f9_min$V2==46113)]<- 46102

f10_min <- fread("flaps12_min_0010.txt")
f10_min$anim <- f10_min$V5 + f10_min$V6
f10_min$V2[which(f10_min$V2==46113)]<- 46102

# Max
f1_max <- fread("flaps12_max_0001.txt") 
f1_max$anim <- f1_max$V5 + f1_max$V6
f1_max$V2[which(f1_max$V2==46113)]<- 46102 

f2_max <- fread("flaps12_max_0002.txt")
f2_max$anim <- f2_max$V5 + f2_max$V6
f2_max$V2[which(f2_max$V2==46113)]<- 46102

f3_max <- fread("flaps12_max_0003.txt")
f3_max$anim <- f3_max$V5 + f3_max$V6
f3_max$V2[which(f3_max$V2==46113)]<- 46102

f4_max <- fread("flaps12_max_0004.txt")
f4_max$anim <- f4_max$V5 + f4_max$V6
f4_max$V2[which(f4_max$V2==46113)]<- 46102

f5_max <- fread("flaps12_max_0005.txt")
f5_max$anim <- f5_max$V5 + f5_max$V6
f5_max$V2[which(f5_max$V2==46113)]<- 46102

f6_max <- fread("flaps12_max_0006.txt")
f6_max$anim <- f6_max$V5 + f6_max$V6
f6_max$V2[which(f6_max$V2==46113)]<- 46102

f7_max <- fread("flaps12_max_0007.txt")
f7_max$anim <- f7_max$V5 + f7_max$V6
f7_max$V2[which(f7_max$V2==46113)]<- 46102

f8_max <- fread("flaps12_max_0008.txt")
f8_max$anim <- f8_max$V5 + f8_max$V6
f8_max$V2[which(f8_max$V2==46113)]<- 46102

f9_max <- fread("flaps12_max_0009.txt")
f9_max$anim <- f9_max$V5 + f9_max$V6
f9_max$V2[which(f9_max$V2==46113)]<- 46102

f10_max <- fread("flaps12_max_0010.txt")
f10_max$anim <- f10_max$V5 + f10_max$V6
f10_max$V2[which(f10_max$V2==46113)]<- 46102

# Number of Animals Infected

# Start the county.anim file

# read the summary file
sum.file <- fread(summary.fnames.90[1], header = TRUE, select = c("Rep","Num_Inf", "Seed_Farms", "Seed_FIPS"))
# Select the right FLAPS file
FLAP <- unlist(strsplit(unlist(strsplit(summary.fnames.90[1],"12_"))[2],"_"))[2] 
flaps <- get(paste0("f",as.numeric(FLAP), "_min"))

# merge FLAPS with summary file to get # animals on seed farm
sum.flaps <- merge(sum.file,flaps[,c("V1", "anim")],by.x="Seed_Farms",by.y="V1",all.x=T)

# Read detail file
det.file <- fread(detail.fnames.90[1], header = TRUE, select = c("Rep","ExposedID", "ExposedCounty", "SourceID", "SourceCounty", "ControlPrevented"))

# Filter on not prevented and unique ExposureID by Rep
det.file <- det.file[det.file$ControlPrevented == "none",]
det.file <- unique(as.data.table(det.file),by=c("Rep","ExposedID")) 

# Summarize the number of animals infected per Rep
rep.totals <- merge(flaps[,c("V1", "V2", "anim")], det.file[,c("Rep","ExposedID", "ExposedCounty")], 
                    by.x = c("V1", "V2"), by.y = c("ExposedID", "ExposedCounty"), all = TRUE)
rep.totals <- as.data.frame(aggregate(rep.totals$anim, by = list(rep.totals$Rep), FUN = sum))
colnames(rep.totals) <- c("Rep", "anim")

# If a Rep's Num_Inf >1, get the additional (non seed farm) animals infected in that rep from the detail file
for(row in 1:max(sum.flaps$Rep)){
  if(sum.flaps$Num_Inf[row]>1){sum.flaps$anim[row]=sum.flaps$anim[row]+rep.totals$anim[which(rep.totals$Rep==sum.flaps$Rep[row])]}
}

county.anim=sum.flaps[,c("Seed_FIPS","anim")]

# Get type
county.anim$Type <- unlist(strsplit(summary.fnames.90[1], "_flaps"))[1]
county.anim$Type <- unlist(strsplit(county.anim$Type, "/"))[2]
county.anim<-as.data.frame(county.anim)



# Read through each summary/detail file 
for (file in 2:length(detail.fnames.90)){
  # read summary file
  sum.file <- fread(summary.fnames.90[file], header = TRUE, select = c("Rep","Num_Inf", "Seed_Farms", "Seed_FIPS"))
  
  # Select the right FLAPS file
  FLAP=unlist(strsplit(unlist(strsplit(summary.fnames.90[file],"12_"))[2],"_"))[2]  
  flaps <- get(paste0("f",as.numeric(FLAP), "_min"))
  
  # merge FLAPS with summary file 
  sum.flaps=merge(sum.file,flaps[,c("V1", "anim")],by.x="Seed_Farms",by.y="V1",all.x=T)
  
  # Read detail file
  det.file <- fread(detail.fnames.90[file], header = TRUE, select = c("Rep","ExposedID", "ExposedCounty", "SourceID", "SourceCounty", "ControlPrevented"))
  
  # Filter on not prevented and unique ExposureID by Rep
  det.file <- det.file[det.file$ControlPrevented == "none",]
  det.file <- unique(as.data.table(det.file),by=c("Rep","ExposedID")) 
  
  # Summarize the number of animals infected per Rep
  rep.totals <- merge(flaps[,c("V1", "V2", "anim")], det.file[,c("Rep","ExposedID", "ExposedCounty")], 
                      by.x = c("V1", "V2"), by.y = c("ExposedID", "ExposedCounty"), all = TRUE)
  rep.totals <- as.data.frame(aggregate(rep.totals$anim, by = list(rep.totals$Rep), FUN = sum))
  colnames(rep.totals) <- c("Rep", "anim")
  
  # If a Rep's Num_Inf >1, get the additional (non seed farm) animals infected in that rep from the detail file
  for(row in 1:max(sum.flaps$Rep)){
    if(sum.flaps$Num_Inf[row]>1){sum.flaps$anim[row]=sum.flaps$anim[row]+rep.totals$anim[which(rep.totals$Rep==sum.flaps$Rep[row])]}
  }
  
  rep.anim=sum.flaps[,c("Seed_FIPS","anim")]
  
  # Get type
  rep.anim$Type <- unlist(strsplit(summary.fnames.90[file], "_flaps"))[1]
  rep.anim$Type <- unlist(strsplit(rep.anim$Type, "/"))[2]
  rep.anim <- as.data.frame(rep.anim)
  
  county.anim=merge(county.anim,rep.anim, by="Seed_FIPS")
  
} # End summary/detail file loop for animals infected

names(county.anim)[seq(3,ncol(county.anim),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(county.anim))/2))))
names(county.anim)[seq(2,ncol(county.anim),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(county.anim))/2))))
#county.anim[,2] <- as.numeric(county.anim[,2])
write.csv(county.anim, "county_anim.csv", row.names = F)

# Convert to long format for plotting
long.county.anim <- reshape(county.anim, idvar = c("Seed_FIPS"), direction = "long", v.names = c("type", "Value"), 
                           varying = list(c(grep("type_", colnames(county.anim))), c(grep("run_", colnames(county.anim)))))
rownames(long.county.anim) <- seq(1:nrow(long.county.anim))
long.county.anim <- long.county.anim[,-2]
long.county.anim <- long.county.anim[!is.na(long.county.anim$Value),]

# Rename the types of runs
long.county.anim$type <- as.factor(long.county.anim$type)
levels(long.county.anim$type)[levels(long.county.anim$type) == "base"] <- "Base"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_MvmtBan_90"] <- "IP Cull, 90% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_MvmtBan_75"] <- "IP Cull, 75% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_DC_MvmtBan_90"] <- "IP & DC Cull, 90% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_DC_MvmtBan_75"] <- "IP & DC Cull, 75% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_MvmtBan_90"] <- "IP Cull & DC Vax, 90% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_MvmtBan_75"] <- "IP Cull & DC Vax, 75% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_3km_MvmtBan_90"] <- "IP Cull & 3km Ring Vax, 90% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_3km_MvmtBan_75"] <- "IP Cull & 3km Ring Vax, 75% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_10km_MvmtBan_90"] <- "IP Cull & 10km Ring Vax, 90% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "IP_VAX_10km_MvmtBan_75"] <- "IP Cull & 10km Ring Vax, 75% Ban"
levels(long.county.anim$type)[levels(long.county.anim$type) == "base_ShipmentsOff"] <- "No Shipments"

color_paired <- c("#000000", "#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                  "#CAB2D6", "#6A3D9A")

head(long.county.anim)
summary(long.county.anim)

# Divide animal numbers by 1000 to print in thousands instead of scientific notation
long.county.anim$Value <- long.county.anim$Value/1000

# Divide the data
long.county.anim.less10 <- long.county.anim[long.county.anim$Value > 10,]
long.county.anim.low <- long.county.anim[long.county.anim$Value <= 10000,]
long.county.anim.high <- long.county.anim[long.county.anim$Value > 10000,]

table(long.county.anim$Value)

jpeg("Anim_Less10_violin_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.anim.less10, aes(x = long.county.anim.less10$type, y = long.county.anim.less10$Value, color = long.county.anim.less10$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") + 
  labs(x = NULL, y = "Animals Infected (thousands)", title = expression("# Infected Animals" > "10 Thousand Subset")) + scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

jpeg("Anim_Low_violin_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.anim.low, aes(x = long.county.anim.low$type, y = long.county.anim.low$Value, color = long.county.anim.low$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + 
  scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Animals Infected (thousands)", title = expression("# Infected Animals" <= "10000 Thousand Subset"))
dev.off()

jpeg("Anim_High_violin_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.anim.high, aes(x = long.county.anim.high$type, y = long.county.anim.high$Value, color = long.county.anim.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")  + ylim(c(10000, 32000)) + 
  scale_color_manual(values = color_paired) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + labs(x = NULL, y = "Animals Infected (thousands)", title = expression("# Infected Animals" > "10000 Thousand Animals"))
dev.off()

jpeg("Anim_High_violin_min_swap.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_violin(data = long.county.anim.high, aes(x = long.county.anim.high$type, y = long.county.anim.high$Value, color = long.county.anim.high$type)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none")+ ylim(c(10000, 32000)) + 
  coord_flip() + scale_color_manual(values = color_paired)  + labs(x = NULL, y = "Animals Infected (thousands)", title = expression("# Infected Animals" > "10000 Thousand Animals"))
dev.off()

## Create individual data frames for each run type
df <- county.anim[,c(1,2:201)]
df1 <- df[df$type_001 != 0 & !is.na(df$type_001),]
newname <- df1$type_001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,202:401)]
df1 <- df[df$type_101 != 0 & !is.na(df$type_101),]
newname <- df1$type_101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,402:601)]
df1 <- df[df$type_201 != 0 & !is.na(df$type_201),]
newname <- df1$type_201[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,602:801)]
df1 <- df[df$type_301 != 0 & !is.na(df$type_301),]
newname <- df1$type_301[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,802:1001)]
df1 <- df[df$type_401 != 0 & !is.na(df$type_401),]
newname <- df1$type_401[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,1002:1201)]
df1 <- df[df$type_501 != 0 & !is.na(df$type_501),]
newname <- df1$type_501[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,1202:1401)]
df1 <- df[df$type_601 != 0 & !is.na(df$type_601),]
newname <- df1$type_601[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,1402:1601)]
df1 <- df[df$type_701 != 0 & !is.na(df$type_701),]
newname <- df1$type_701[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,1602:1801)]
df1 <- df[df$type_801 != 0 & !is.na(df$type_801),]
newname <- df1$type_801[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,1802:2001)]
df1 <- df[df$type_901 != 0 & !is.na(df$type_901),]
newname <- df1$type_901[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,2002:2201)]
df1 <- df[df$type_1001 != 0 & !is.na(df$type_1001),]
newname <- df1$type_1001[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)

df <- county.anim[,c(1,2202:2401)]
df1 <- df[df$type_1101 != 0 & !is.na(df$type_1101),]
newname <- df1$type_1101[1]
df.new <- df[,-grep("type_", colnames(df))]
assign(paste0(newname, "_Anim"), df.new)


## calculate median of each run
base_Anim$median <- apply(base_Anim[grep("run_", colnames(base_Anim))], 1, median, na.rm = TRUE)/1000
IP_MvmtBan_90_Anim$median <- apply(IP_MvmtBan_90_Anim[grep("run_", colnames(IP_MvmtBan_90_Anim))], 1, median, na.rm = TRUE)/1000
IP_MvmtBan_75_Anim$median <- apply(IP_MvmtBan_75_Anim[grep("run_", colnames(IP_MvmtBan_75_Anim))], 1, median, na.rm = TRUE)/1000
IP_DC_MvmtBan_90_Anim$median <- apply(IP_DC_MvmtBan_90_Anim[grep("run_", colnames(IP_DC_MvmtBan_90_Anim))], 1, median, na.rm = TRUE)/1000
IP_DC_MvmtBan_75_Anim$median <- apply(IP_DC_MvmtBan_75_Anim[grep("run_", colnames(IP_DC_MvmtBan_75_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_MvmtBan_90_Anim$median <- apply(IP_VAX_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_MvmtBan_90_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_MvmtBan_75_Anim$median <- apply(IP_VAX_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_MvmtBan_75_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_3km_MvmtBan_90_Anim$median <- apply(IP_VAX_3km_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_3km_MvmtBan_75_Anim$median <- apply(IP_VAX_3km_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_10km_MvmtBan_90_Anim$median <- apply(IP_VAX_10km_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_Anim))], 1, median, na.rm = TRUE)/1000
IP_VAX_10km_MvmtBan_75_Anim$median <- apply(IP_VAX_10km_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_Anim))], 1, median, na.rm = TRUE)/1000
base_ShipmentsOff_Anim$median <- apply(base_ShipmentsOff_Anim[grep("run_", colnames(base_ShipmentsOff_Anim))], 1, median, na.rm = TRUE)

# Create a data.frame with just the fips code and median to map
base_Anim_med <- base_Anim[,c(1,102)]
IP_MvmtBan_90_Anim_med <- IP_MvmtBan_90_Anim[,c(1,102)]
IP_MvmtBan_75_Anim_med <- IP_MvmtBan_75_Anim[,c(1,102)]
IP_DC_MvmtBan_90_Anim_med <- IP_DC_MvmtBan_90_Anim[,c(1,102)]
IP_DC_MvmtBan_75_Anim_med <- IP_DC_MvmtBan_75_Anim[,c(1,102)]
IP_VAX_MvmtBan_90_Anim_med <- IP_VAX_MvmtBan_90_Anim[,c(1,102)]
IP_VAX_MvmtBan_75_Anim_med <- IP_VAX_MvmtBan_75_Anim[,c(1,102)]
IP_VAX_3km_MvmtBan_90_Anim_med <- IP_VAX_3km_MvmtBan_90_Anim[,c(1,102)]
IP_VAX_3km_MvmtBan_75_Anim_med <- IP_VAX_3km_MvmtBan_75_Anim[,c(1,102)]
IP_VAX_10km_MvmtBan_90_Anim_med <- IP_VAX_10km_MvmtBan_90_Anim[,c(1,102)]
IP_VAX_10km_MvmtBan_75_Anim_med <- IP_VAX_10km_MvmtBan_75_Anim[,c(1,102)]
base_ShipmentsOff_Anim_med <- base_ShipmentsOff_Anim[,c(1,102)]


## Calculate the Upper 2.5% for each county
base_Anim$upper <- apply(base_Anim[grep("run_", colnames(base_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_MvmtBan_90_Anim$upper <- apply(IP_MvmtBan_90_Anim[grep("run_", colnames(IP_MvmtBan_90_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_MvmtBan_75_Anim$upper <- apply(IP_MvmtBan_75_Anim[grep("run_", colnames(IP_MvmtBan_75_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_DC_MvmtBan_90_Anim$upper <- apply(IP_DC_MvmtBan_90_Anim[grep("run_", colnames(IP_DC_MvmtBan_90_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_DC_MvmtBan_75_Anim$upper <- apply(IP_DC_MvmtBan_75_Anim[grep("run_", colnames(IP_DC_MvmtBan_75_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_MvmtBan_90_Anim$upper <- apply(IP_VAX_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_MvmtBan_90_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_MvmtBan_75_Anim$upper <- apply(IP_VAX_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_MvmtBan_75_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_3km_MvmtBan_90_Anim$upper <- apply(IP_VAX_3km_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_3km_MvmtBan_90_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_3km_MvmtBan_75_Anim$upper <- apply(IP_VAX_3km_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_3km_MvmtBan_75_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_10km_MvmtBan_90_Anim$upper <- apply(IP_VAX_10km_MvmtBan_90_Anim[grep("run_", colnames(IP_VAX_10km_MvmtBan_90_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
IP_VAX_10km_MvmtBan_75_Anim$upper <- apply(IP_VAX_10km_MvmtBan_75_Anim[grep("run_", colnames(IP_VAX_10km_MvmtBan_75_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000
base_ShipmentsOff_Anim$upper <- apply(base_ShipmentsOff_Anim[grep("run_", colnames(base_ShipmentsOff_Anim))], 1, quantile, probs=0.975, na.rm = TRUE)/1000

# Create a data.frame with just the fips code and upper to map
base_Anim_up <- base_Anim[,c(1,103)]
IP_MvmtBan_90_Anim_up <- IP_MvmtBan_90_Anim[,c(1,103)]
IP_MvmtBan_75_Anim_up <- IP_MvmtBan_75_Anim[,c(1,103)]
IP_DC_MvmtBan_90_Anim_up <- IP_DC_MvmtBan_90_Anim[,c(1,103)]
IP_DC_MvmtBan_75_Anim_up <- IP_DC_MvmtBan_75_Anim[,c(1,103)]
IP_VAX_MvmtBan_90_Anim_up <- IP_VAX_MvmtBan_90_Anim[,c(1,103)]
IP_VAX_MvmtBan_75_Anim_up <- IP_VAX_MvmtBan_75_Anim[,c(1,103)]
IP_VAX_3km_MvmtBan_90_Anim_up <- IP_VAX_3km_MvmtBan_90_Anim[,c(1,103)]
IP_VAX_3km_MvmtBan_75_Anim_up <- IP_VAX_3km_MvmtBan_75_Anim[,c(1,103)]
IP_VAX_10km_MvmtBan_90_Anim_up <- IP_VAX_10km_MvmtBan_90_Anim[,c(1,103)]
IP_VAX_10km_MvmtBan_75_Anim_up <- IP_VAX_10km_MvmtBan_75_Anim[,c(1,103)]
base_ShipmentsOff_Anim_up <- base_ShipmentsOff_Anim[,c(1,103)]

# Create a vector for scale
Anim_scale <- c(base_Anim$upper, IP_MvmtBan_90_Anim$upper, IP_MvmtBan_75_Anim$upper, IP_DC_MvmtBan_90_Anim$upper, 
               IP_DC_MvmtBan_75_Anim$upper, IP_VAX_MvmtBan_90_Anim$upper, IP_VAX_MvmtBan_75_Anim$upper, 
               IP_VAX_3km_MvmtBan_90_Anim$upper, IP_VAX_3km_MvmtBan_75_Anim$upper, IP_VAX_10km_MvmtBan_90_Anim$upper, 
               IP_VAX_10km_MvmtBan_75_Anim$upper, base_ShipmentsOff_Anim$upper)
Anim_scale <- round(quantile(Anim_scale, probs = c(0, .1, .2, .3, .4, .5, .6, .7, 0.75, .8, 0.85, .9, 0.95, 0.99, 1), na.rm = TRUE), 0)
Anim_scale <- c(0, 1, 2, 10, 100, 1000, 5000, 10000, 20000, 25000, 30000, 32000)

# Median maps
jpeg("Anim_Med_Base_max.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Med_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_Anim_med, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

# Upper 2.5% maps
jpeg("Anim_Up_Base_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IP90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_90_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IP75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_MvmtBan_75_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPDC90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_90_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPDC75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_DC_MvmtBan_75_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX90_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_90_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX75_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_MvmtBan_75_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX90_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_90_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX75_3_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_3km_MvmtBan_75_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX90_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_90_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_IPVAX75_10_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(IP_VAX_10km_MvmtBan_75_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

jpeg("Anim_Up_ShipmentsOff_min.jpeg", width = 760, height = 520, units = 'px', res = 100)
map_by_fips(base_ShipmentsOff_Anim_up, county.border.col = NA, state.border.col = "gray30", missing.include = TRUE, color.break.type = "values", 
            color.break.values = Anim_scale, color.sequence = color_YlOrRd, legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
dev.off()

#### Routes of Infection (Detail Files) ####

detail.fname.90 <- detail.fnames.90[1]
county.rt <- as.data.frame(fread(detail.fname.90, header = TRUE, select = c("SourceCounty", "ExposedCounty", "ControlPrevented")))
detail.res <- detail.res[detail.res$ControlPrevented == "none",]
detail.res <- detail.res[,c("SourceCounty", "ExposedCounty")]
county.rt$Type <- unlist(strsplit(detail.fname.90, "_flaps"))[1]
county.rt$Type <- unlist(strsplit(county.rt$Type, "/"))[2]
county.rt <- county.rt[county.rt$SourceCounty != county.rt$ExposedCounty,]
county.rt$Freq <- rep(1, nrow(county.rt))
county.rt <- aggregate(county.rt$Freq, by = list(SourceCounty = county.rt$SourceCounty, ExposedCounty = county.rt$ExposedCounty, Type = county.rt$Type), FUN = sum)
colnames(county.rt)[4] <- "Freq"


for(i in 2:length(detail.fnames.90)){
  detail.fname.90 <- detail.fnames.90[i]
  detail.res <- as.data.frame(fread(detail.fname.90, header = TRUE, select = c("SourceCounty", "ExposedCounty", "ControlPrevented")))
  detail.res <- detail.res[detail.res$ControlPrevented == "none",]
  detail.res <- detail.res[,c("SourceCounty", "ExposedCounty")]
  detail.res$Type <- unlist(strsplit(detail.fname.90, "_flaps"))[1]
  detail.res$Type <- unlist(strsplit(detail.res$Type, "/"))[2]
  detail.res <- detail.res[detail.res$SourceCounty != detail.res$ExposedCounty,]
  detail.res$Freq <- rep(1, nrow(detail.res))
  detail.res <- as.data.frame(detail.res)
  county.rt <- rbind(county.rt, detail.res)
  county.rt <- aggregate(county.rt$Freq, by = list(SourceCounty = county.rt$SourceCounty, ExposedCounty = county.rt$ExposedCounty, Type = county.rt$Type), FUN = sum)
  colnames(county.rt)[4] <- "Freq"
}

write.csv(county.rt, "county_rt.csv")
# county.rt <- read.table("county_rt.csv", header = TRUE)

## Create individual data frames for each run type
base_rt <- county.rt[county.rt$Type == "base_",]
IP_MvmtBan_90_rt <- county.rt[county.rt$Type == "IP_MvmtBan_90_",]
IP_DC_MvmtBan_90_rt <- county.rt[county.rt$Type == "IP_DC_MvmtBan_90_",]
base_ShipmentsOff_rt <- county.rt[county.rt$Type == "base_ShipmentsOff",]
IP_MvmtBan_75_rt <- county.rt[county.rt$Type == "IP_MvmtBan_75",]
IP_DC_MvmtBan_75_rt <- county.rt[county.rt$Type == "IP_DC_MvmtBan_75",]

# Merge the base and each control run
IP_90_diff_rts <- merge(IP_MvmtBan_90_rt, base_rt, by = c("SourceCounty", "NewCounty"), all = TRUE) 
IP_DC_90_diff_rts <- merge(IP_DC_MvmtBan_90_rt, base_rt, by = c("SourceCounty", "NewCounty"), all = TRUE) 
IPVAX_90_diff_rts <- merge(base_ShipmentsOff_rt, base_rt, by = c("SourceCounty", "NewCounty"), all = TRUE)
IPVAX_90_3_diff_rts <- merge(IP_MvmtBan_75_rt, base_rt, by = c("SourceCounty", "NewCounty"), all = TRUE)
IPVAX_90_10_diff_rts <- merge(IP_DC_MvmtBan_75_rt, base_rt, by = c("SourceCounty", "NewCounty"), all = TRUE)

# Reduce the number of columns
IP_90_diff_rts <- IP_90_diff_rts[,c(1,2,4,6)]
IP_DC_90_diff_rts <- IP_DC_90_diff_rts[,c(1,2,4,6)]
IPVAX_90_diff_rts <- IPVAX_90_diff_rts[,c(1,2,4,6)]
IPVAX_90_3_diff_rts <- IPVAX_90_3_diff_rts[,c(1,2,4,6)]
IPVAX_90_10_diff_rts <- IPVAX_90_10_diff_rts[,c(1,2,4,6)]

# Create a column for the difference between base and control
colnames(IP_90_diff_rts) <- c("SourceCounty", "NewCounty", "CtrlFreq", "BaseFreq")
IP_90_diff_rts$diff <- IP_90_diff_rts$BaseFreq - IP_90_diff_rts$CtrlFreq
IP_90_diff_rts <- IP_90_diff_rts[which(!is.na(IP_90_diff_rts$diff)),]
colnames(IP_DC_90_diff_rts) <- c("SourceCounty", "NewCounty", "CtrlFreq", "BaseFreq")
IP_DC_90_diff_rts$diff <- IP_DC_90_diff_rts$BaseFreq - IP_DC_90_diff_rts$CtrlFreq
IP_DC_90_diff_rts <- IP_DC_90_diff_rts[which(!is.na(IP_DC_90_diff_rts$diff)),]
colnames(IPVAX_90_diff_rts) <- c("SourceCounty", "NewCounty", "CtrlFreq", "BaseFreq")
IPVAX_90_diff_rts$diff <- IPVAX_90_diff_rts$BaseFreq - IPVAX_90_diff_rts$CtrlFreq
IPVAX_90_diff_rts <- IPVAX_90_diff_rts[which(!is.na(IPVAX_90_diff_rts$diff)),]
colnames(IPVAX_90_3_diff_rts) <- c("SourceCounty", "NewCounty", "CtrlFreq", "BaseFreq")
IPVAX_90_3_diff_rts$diff <- IPVAX_90_3_diff_rts$BaseFreq - IPVAX_90_3_diff_rts$CtrlFreq
IPVAX_90_3_diff_rts <- IPVAX_90_3_diff_rts[which(!is.na(IPVAX_90_3_diff_rts$diff)),]
colnames(IPVAX_90_10_diff_rts) <- c("SourceCounty", "NewCounty", "CtrlFreq", "BaseFreq")
IPVAX_90_10_diff_rts$diff <- IPVAX_90_10_diff_rts$BaseFreq - IPVAX_90_10_diff_rts$CtrlFreq
IPVAX_90_10_diff_rts <- IPVAX_90_10_diff_rts[which(!is.na(IPVAX_90_10_diff_rts$diff)),]

### Most Controlled Routes of Infection ####

Ctrlled_IP_90 <- head(arrange(IP_90_diff_rts, desc(diff)), n = 500)
Ctrlled_IP_DC_90 <- head(arrange(IP_DC_90_diff_rts, desc(diff)), n = 500)
Ctrlled_IPVAX_90 <- head(arrange(IPVAX_90_diff_rts, desc(diff)), n = 500)
Ctrlled_IPVAX_90_3 <- head(arrange(IPVAX_90_3_diff_rts, desc(diff)), n = 500)
Ctrlled_IPVAX_90_10 <- head(arrange(IPVAX_90_10_diff_rts, desc(diff)), n = 500)

# Get lat/long data
county <- map_data("county")
state <-  map_data("state")
county_data <- read.csv("county_data.csv", header = TRUE)
Ctrlled_IP_90 <- merge(Ctrlled_IP_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(Ctrlled_IP_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
Ctrlled_IP_90 <- merge(Ctrlled_IP_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(Ctrlled_IP_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

Ctrlled_IP_DC_90 <- merge(Ctrlled_IP_DC_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(Ctrlled_IP_DC_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
Ctrlled_IP_DC_90 <- merge(Ctrlled_IP_DC_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(Ctrlled_IP_DC_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

Ctrlled_IPVAX_90 <- merge(Ctrlled_IPVAX_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
Ctrlled_IPVAX_90 <- merge(Ctrlled_IPVAX_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

Ctrlled_IPVAX_90_3 <- merge(Ctrlled_IPVAX_90_3, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90_3) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
Ctrlled_IPVAX_90_3 <- merge(Ctrlled_IPVAX_90_3, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90_3) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

Ctrlled_IPVAX_90_10 <- merge(Ctrlled_IPVAX_90_10, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90_10) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
Ctrlled_IPVAX_90_10 <- merge(Ctrlled_IPVAX_90_10, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(Ctrlled_IPVAX_90_10) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

min_ctrlled <- min(c(Ctrlled_IP_90$Diff, Ctrlled_IP_DC_90$Diff, Ctrlled_IPVAX_90$Diff, Ctrlled_IPVAX_90_3$Diff, Ctrlled_IPVAX_90_10$Diff))
max_ctrlled <- max(c(Ctrlled_IP_90$Diff, Ctrlled_IP_DC_90$Diff, Ctrlled_IPVAX_90$Diff, Ctrlled_IPVAX_90_3$Diff, Ctrlled_IPVAX_90_10$Diff))


# Maps
jpeg("Rt_IP90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = Ctrlled_IP_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = rev(as.numeric(Diff)), 
                            size = 1, curvature = 0.1)) + scale_color_gradientn(colours = color_purple[c(5:8)], limits = c(min_ctrlled, max_ctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("Rt_IP_DC90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = Ctrlled_IP_DC_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = color_purple[c(5:8)], limits = c(min_ctrlled, max_ctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("Rt_IPVAX90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = Ctrlled_IPVAX_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = color_purple[c(5:8)], limits = c(min_ctrlled, max_ctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("Rt_IPVAX90_3.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = Ctrlled_IPVAX_90_3, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = color_purple[c(5:8)], limits = c(min_ctrlled, max_ctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("Rt_IPVAX90_10.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = Ctrlled_IPVAX_90_3, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = color_purple[c(5:8)], limits = c(min_ctrlled, max_ctrlled)) + 
  theme(legend.position = "none")
dev.off()


### Least Controlled Routes of Infection ####
NotCtrlled_IP_90 <- head(arrange(IP_90_diff_rts, diff), n = 500)
NotCtrlled_IP_DC_90 <- head(arrange(IP_DC_90_diff_rts, diff), n = 500)
NotCtrlled_IPVAX_90 <- head(arrange(IPVAX_90_diff_rts, diff), n = 500)
NotCtrlled_IPVAX_90_3 <- head(arrange(IPVAX_90_3_diff_rts, diff), n = 500)
NotCtrlled_IPVAX_90_10 <- head(arrange(IPVAX_90_10_diff_rts, diff), n = 500)

# Get lat/long data
county <- map_data("county")
state <-  map_data("state")
county_data <- read.csv("county_data.csv", header = TRUE)
NotCtrlled_IP_90 <- merge(NotCtrlled_IP_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(NotCtrlled_IP_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
NotCtrlled_IP_90 <- merge(NotCtrlled_IP_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(NotCtrlled_IP_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

NotCtrlled_IP_DC_90 <- merge(NotCtrlled_IP_DC_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(NotCtrlled_IP_DC_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
NotCtrlled_IP_DC_90 <- merge(NotCtrlled_IP_DC_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(NotCtrlled_IP_DC_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

NotCtrlled_IPVAX_90 <- merge(NotCtrlled_IPVAX_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
NotCtrlled_IPVAX_90 <- merge(NotCtrlled_IPVAX_90, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

NotCtrlled_IPVAX_90_3 <- merge(NotCtrlled_IPVAX_90_3, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90_3) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
NotCtrlled_IPVAX_90_3 <- merge(NotCtrlled_IPVAX_90_3, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90_3) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

NotCtrlled_IPVAX_90_10 <- merge(NotCtrlled_IPVAX_90_10, county_data[,c("FIPS", "LONG", "LAT")], by.x = "SourceCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90_10) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat")
NotCtrlled_IPVAX_90_10 <- merge(NotCtrlled_IPVAX_90_10, county_data[,c("FIPS", "LONG", "LAT")], by.x = "NewCounty", by.y = "FIPS")
colnames(NotCtrlled_IPVAX_90_10) <- c("SourceCounty", "NewCounty", "CFreq", "BFreq", "Diff", "S.long", "S.lat", "N.long", "N.lat")

min_notctrlled <- min(c(NotCtrlled_IP_90$Diff, NotCtrlled_IP_DC_90$Diff, NotCtrlled_IPVAX_90$Diff, NotCtrlled_IPVAX_90_3$Diff, NotCtrlled_IPVAX_90_10$Diff))
max_notctrlled <- max(c(NotCtrlled_IP_90$Diff, NotCtrlled_IP_DC_90$Diff, NotCtrlled_IPVAX_90$Diff, NotCtrlled_IPVAX_90_3$Diff, NotCtrlled_IPVAX_90_10$Diff))


# Maps
jpeg("RtNC_IP90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = NotCtrlled_IP_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = rev(color_orange[c(5:8)]), limits = c(min_notctrlled, max_notctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("RtNC_IP_DC90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = NotCtrlled_IP_DC_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = rev(color_orange[c(5:8)]), limits = c(min_notctrlled, max_notctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("RtNC_IPVAX90.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = NotCtrlled_IPVAX_90, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = rev(color_orange[c(5:8)]), limits = c(min_notctrlled, max_notctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("RtNC_IPVAX90_3.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = NotCtrlled_IPVAX_90_3, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = rev(color_orange[c(5:8)]), limits = c(min_notctrlled, max_notctrlled)) + 
  theme(legend.position = "none")
dev.off()

jpeg("RtNC_IPVAX90_10.jpeg", width = 760, height = 520, units = 'px', res = 100)
ggplot() + geom_polygon(data = county, aes(x = county$long, y = county$lat, group = county$group), color = "gray", fill = NA) + 
  coord_fixed(1.3) + geom_polygon(data = state, aes(x = long, y = lat, group = group), color = "black", fill = NA) + 
  theme_void() + geom_curve(data = NotCtrlled_IPVAX_90_3, aes(x = S.long, y = S.lat, xend = N.long, yend = N.lat, color = as.numeric(Diff)), 
                            size = 1, curvature = 0.1) + scale_color_gradientn(colours = rev(color_orange[c(5:8)]), limits = c(min_notctrlled, max_notctrlled)) + 
  labs(title = "Bottom IPVAX10") + theme(legend.position = "none")
dev.off()

#### Create a Summary Table ####

Summary_90 <- cbind(NI_Medians_90_Low, NI_Medians_90_High, IP_Medians_90_Low, IP_Medians_90_High, DC_Medians_90_Low, DC_Medians_90_High,
                    Vacc_Medians_90_Low, Vacc_Medians_90_High, Dur_Medians_90_Low, Dur_Medians_90_High, round(Prop_90, 1))
colnames(Summary_90) <- c("</= 1000", "> 1000", "</= 3000", "> 3000", "</= 2000", "> 2000", 
                          "</= 3000", "> 3000", "</= 100", "> 100", "> 10 Infections")
rownames(Summary_90) <- c("Base", "IP Culls", "IP & DC Culls", "IP Culls & DC Vaccination", "IP Culls & 3 km Vaccination Ring", 
                          "IP Culls & 10 km Vaccination Ring")

jpeg("Summary_90.jpeg", width = 793, height = 400, units = 'px', res = 100)
kable(Summary_90, "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), full_width = F, position = "float_right")%>%
  add_header_above(c(" ", "Number of IPs" = 2, "IPs Culled" = 2, "DCs Culled" = 2, "Vaccine Doses" = 2, 
                     "Duration (days)" = 2, "% of Runs" = 1))%>%
  add_header_above(c(" ", "90% Shipment Ban - Median Values" = 10, ""))


dev.off()

Summary_75 <- cbind(NI_Medians_75_Low, NI_Medians_75_High, IP_Medians_75_Low, IP_Medians_75_High, DC_Medians_75_Low, DC_Medians_75_High,
                    Vacc_Medians_75_Low, Vacc_Medians_75_High, Dur_Medians_75_Low, Dur_Medians_75_High, round(Prop_75, 1))
colnames(Summary_75) <- c("</= 1000", "> 1000", "</= 3000", "> 3000", "</= 2000", "> 2000", 
                          "</= 3000", "> 3000", "</= 100", "> 100", "> 10 Infections")
rownames(Summary_75) <- c("Base", "IP Culls", "IP & DC Culls", "IP Culls & DC Vaccination", "IP Culls & 3 km Vaccination Ring", 
                          "IP Culls & 10 km Vaccination Ring")

jpeg("Summary_75.jpeg", width = 793, height = 400, units = 'px', res = 100)
kable(Summary_75, "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), full_width = F, position = "float_right")%>%
  add_header_above(c(" ", "Number of IPs" = 2, "IPs Culled" = 2, "DCs Culled" = 2, "Vaccine Doses" = 2, 
                     "Duration (days)" = 2, "% of Runs" = 1))%>%
  add_header_above(c(" ", "75% Shipment Ban - Median Values" = 10, ""))


dev.off()



#### Combined Violin Plot

jpeg("NI_Combined.jpeg", width = 760, height = 520, units = 'px', res = 100)
dev.off()

low <- ggplot() + geom_violin(data = long.county.inf.low, aes(x = long.county.inf.low$type, y = long.county.inf.low$Value)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_bw() + theme(legend.position = "none") +
  labs(x = NULL, y = NULL,  title = NULL) + 
  theme(axis.text.x = element_text(size = 14)) + coord_flip() + ylim(c(1, 5000)) + annotate("text", x = 3000, y = 2, label = "3.58 million simulations") + 
  geom_dotplot(data = long.county.inf.low, aes(x = long.county.inf.low$type, y = long.county.inf.low$Value, color = long.county.inf.low$Value)) + 
                 scale_color_manual(values = color_YlOrRd)  


high <- ggplot() + geom_violin(data = long.county.inf.high, aes(x = long.county.inf.high$type, y = long.county.inf.high$Value)) + 
  scale_x_discrete(limits = c("Base", "IP Cull, 90% Ban", "IP Cull, 75% Ban", "IP & DC Cull, 90% Ban", "IP & DC Cull, 75% Ban", "IP Cull & DC Vax, 90% Ban", "IP Cull & DC Vax, 75% Ban","IP Cull & 3km Ring Vax, 90% Ban", 
                              "IP Cull & 3km Ring Vax, 75% Ban", "IP Cull & 10km Ring Vax, 90% Ban", "IP Cull & 10km Ring Vax, 75% Ban", "No Shipments")) + theme_minimal() + theme(legend.position = "none") +
  labs(x = NULL, y = NULL,  title = NULL) + annotate("text", x = 5000, y = 2, label = "74 Thousand simulations") +
  theme(axis.text.x = element_blank(), axis.line.x = element_line(color="gray40", linetype = "dashed")) + ylim(c(5000, 32000)) + coord_flip() + 
  geom_dotplot(data = long.county.inf.high, aes(x = long.county.inf.high$type, y = long.county.inf.high$Value, color = long.county.inf.high$Value)) +
  scale_color_manual(values = color_YlOrRd)  

grid.arrange(low, high, ncol=2, top = "Number of Premises Infected")

