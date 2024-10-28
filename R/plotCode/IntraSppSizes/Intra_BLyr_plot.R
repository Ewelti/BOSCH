##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(scales)
library(stringr)
library(rphylopic)
library(vegan)

# attach data
intra <- read.csv("RawData/IntraSppBS.csv", header=T)
head(intra)
nrow(intra)
unique(intra$Species)

#fix some problem where these variables are not numeric
intra$BL <- as.numeric(intra$Body_Length) 
intra$HW <- as.numeric(intra$Head_Width)
intra$BW <- as.numeric(intra$Body_Width)
intra$BH <- as.numeric(intra$Height)
intra$LA <- as.numeric(intra$Length.of.1st.Antennae)

##scale everything
#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
intra <- scaleVars(intra)
head(intra)

#log densities
intra$dens <- as.numeric(intra$density_per_m2)
intra$Ldens <- log10(intra$dens+1)

#get initial summary data
library(dplyr)
# Summarize the data by 'Species'
summary_intra <- intra %>%
  group_by(Species) %>%
  summarise(
    count = n(),                       # Count of rows per species
    mean_body_length = mean(Body_Length, na.rm = TRUE),     # Mean body length per species
    mean_head_width = mean(Head_Width, na.rm = TRUE),) %>%  # Mean head width per species
  arrange(mean_body_length)
detach("package:dplyr", unload = TRUE)
summary_intra

#calculate BL estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_s <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "BL"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(BL ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_s <- rbind(ests_s, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

ests_s[c('spp','site', 'yr')] <- str_split_fixed(ests_s$spp_site_yr, ' ', 3)
colnames(ests_s)[2] ="BL_est"
colnames(ests_s)[3] ="BL_SE"
ests_s$yr <-as.numeric(ests_s$yr)
ests_s = subset(ests_s, select = -c(spp_site_yr))
head(ests_s)
##

#subset by spp
ed <- ests_s[which(ests_s$spp=="ED"), ]
po <- ests_s[which(ests_s$spp=="PO"), ]
hs <- ests_s[which(ests_s$spp=="HS"), ]
ov <- ests_s[which(ests_s$spp=="OV"), ]
gr <- ests_s[which(ests_s$spp=="GR"), ]
et <- ests_s[which(ests_s$spp=="ET"), ]
af <- ests_s[which(ests_s$spp=="AF"), ]
br <- ests_s[which(ests_s$spp=="BR"), ]
aa <- ests_s[which(ests_s$spp=="AA"), ]


#no sci notation
options(scipen = 999)
options(na.action = "na.omit")

# Function to add silhouette to existing plot
add_phylopic <- function (img = NULL, name = NULL, uuid = NULL, upload_img = NULL, filter = NULL, 
                          x = NULL, y = NULL, ysize = NULL, height = NULL, 
                          width = NULL, alpha = 1, color = NA, fill = "black", 
                          horizontal = FALSE, vertical = FALSE, angle = 0, 
                          hjust = 0.5, vjust = 0.5, remove_background = TRUE, 
                          verbose = FALSE) {
  # Install and load required packages
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    install.packages("rsvg")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    install.packages("grid")
  }
  if (!requireNamespace("lifecycle", quietly = TRUE)) {
    install.packages("lifecycle")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    install.packages("png")
  }
  library(rsvg)
  library(grid)
  library(lifecycle)
  library(png)
  
  # Check for required arguments
  message("Checking required arguments...")
  if (all(sapply(list(img, name, uuid, upload_img), is.null))) {
    stop("One of `img`, `name`, `uuid`, or `upload_img` is required.")
  }
  if (sum(sapply(list(img, name, uuid, upload_img), is.null)) < 3) {
    stop("Only one of `img`, `name`, `uuid`, or `upload_img` may be specified")
  }
  
  # Load image from file if `upload_img` is specified and convert to raster
  if (!is.null(upload_img)) {
    message("Loading image from file...")
    # Convert SVG to PNG and read it as a raster array
    temp_file <- tempfile(fileext = ".png")
    rsvg::rsvg_png(upload_img, file = temp_file)
    img <- png::readPNG(temp_file)
    unlink(temp_file)  # Remove the temporary file after reading
    message("Image successfully loaded and converted.")
  }
  
  # Validation checks
  if (any(alpha > 1 | alpha < 0)) {
    stop("`alpha` must be between 0 and 1.")
  }
  if (any(hjust > 1 | hjust < 0)) {
    stop("`hjust` must be between 0 and 1.")
  }
  if (any(vjust > 1 | vjust < 0)) {
    stop("`vjust` must be between 0 and 1.")
  }
  if (!is.logical(verbose)) {
    stop("`verbose` should be a logical value.")
  }
  
  # Deprecation warning for `ysize`
  if (!is.null(ysize)) {
    lifecycle::deprecate_warn("1.5.0", "add_phylopic(ysize)", 
                              "add_phylopic(height)")
    if (is.null(height)) 
      height <- ysize
  }
  if (!is.null(height) & !is.null(width)) {
    stop("At least one of `height` or `width` must be NULL.")
  }
  
  # Define plot boundaries
  message("Defining plot boundaries...")
  usr <- par()$usr
  usr_x <- if (par()$xlog) 10^usr[1:2] else usr[1:2]
  usr_y <- if (par()$ylog) 10^usr[3:4] else usr[3:4]
  
  # Set default x and y positions if not specified
  if (is.null(x)) {
    mn <- mean(usr[1:2])
    x <- if (par()$xlog) 10^mn else mn
  }
  if (is.null(y)) {
    mn <- mean(usr[3:4])
    y <- if (par()$ylog) 10^mn else mn
  }
  
  # Set default height and width if not specified
  if (is.null(height) && is.null(width)) {
    height <- abs(diff(usr_y)) * 0.1  # Default to 10% of y-axis range if not set
    width <- NA
  } else {
    if (!is.null(height)) {
      base_y <- grconvertY(ifelse(par()$ylog, 1, 0), to = "ndc")
      height <- grconvertY(height, to = "ndc") - base_y
    }
    if (!is.null(width)) {
      base_x <- grconvertX(ifelse(par()$xlog, 1, 0), to = "ndc")
      width <- grconvertX(width, to = "ndc") - base_x
    }
  }
  
  message(paste("Calculated width:", width))
  message(paste("Calculated height:", height))
  
  # Convert x and y coordinates
  message("Converting coordinates...")
  x <- grconvertX(x, to = "ndc")
  y <- grconvertY(y, to = "ndc")
  
  # Render the image
  message("Rendering the image...")
  if (is.null(img)) {
    stop("No image data available for rendering.")
  }
  if (horizontal || vertical) img <- flip_phylopic(img, horizontal, vertical)
  if (angle != 0) img <- rotate_phylopic(img, angle)
  if (is.na(color) || color == "original") color <- NULL
  if (is.na(fill)) {
    fill <- color
    color <- NULL
  }
  if (fill == "original") fill <- NULL
  img <- recolor_phylopic(img, alpha, color, fill, remove_background)
  
  # Display the image as a raster
  grid.raster(img, x = x, y = y, width = width, height = height, 
              hjust = hjust, vjust = vjust)
  message("Image rendering completed.")
}

##BL plot##
tiff(filename = "plots/Year_BL_wSilhouettes.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')
par(mar=c(2.5,5,4,0.5),mfrow=c(3,3))

#####Aphelocheirus aestivalis
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(aa$BL_est-aa$BL_SE-1),max(aa$BL_est + aa$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$BL_est[aa$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$BL_est[aa$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$BL_est[aa$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$BL_est[aa$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$BL_est[aa$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$BL_est[aa$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$BL_est[aa$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$BL_est[aa$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(aa$yr[aa$site=="Auba"]-0.15, aa$BL_est[aa$site=="Auba"]-aa$BL_SE[aa$site=="Auba"], aa$yr[aa$site=="Auba"]-0.15, aa$BL_est[aa$site=="Auba"]+aa$BL_SE[aa$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="Bieb"]-0.05, aa$BL_est[aa$site=="Bieb"]-aa$BL_SE[aa$site=="Bieb"], aa$yr[aa$site=="Bieb"]-0.05, aa$BL_est[aa$site=="Bieb"]+aa$BL_SE[aa$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="O3"]+0.05, aa$BL_est[aa$site=="O3"]-aa$BL_SE[aa$site=="O3"], aa$yr[aa$site=="O3"]+0.05, aa$BL_est[aa$site=="O3"]+aa$BL_SE[aa$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="W1"]+0.15, aa$BL_est[aa$site=="W1"]-aa$BL_SE[aa$site=="W1"], aa$yr[aa$site=="W1"]+0.15, aa$BL_est[aa$site=="W1"]+aa$BL_SE[aa$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/aphelocheirus_aestivalis.svg",
  x = 2002,
  y = min(aa$BL_est - aa$BL_SE - 1) + 1.1,
  width = 2.5,
  height = NULL
)
title("a. Aphelocheirus aestivalis",bty="n",cex.main=2)

###############Ancylus fluviatilis
af[is.na(af)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(af$BL_est-af$BL_SE-1),max(af$BL_est + af$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$BL_est[af$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$BL_est[af$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$BL_est[af$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$BL_est[af$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$BL_est[af$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$BL_est[af$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$BL_est[af$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$BL_est[af$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(af$yr[af$site=="Auba"]-0.15, af$BL_est[af$site=="Auba"]-af$BL_SE[af$site=="Auba"], af$yr[af$site=="Auba"]-0.15, af$BL_est[af$site=="Auba"]+af$BL_SE[af$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="Bieb"]-0.05, af$BL_est[af$site=="Bieb"]-af$BL_SE[af$site=="Bieb"], af$yr[af$site=="Bieb"]-0.05, af$BL_est[af$site=="Bieb"]+af$BL_SE[af$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="O3"]+0.05, af$BL_est[af$site=="O3"]-af$BL_SE[af$site=="O3"], af$yr[af$site=="O3"]+0.05, af$BL_est[af$site=="O3"]+af$BL_SE[af$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="W1"]+0.15, af$BL_est[af$site=="W1"]-af$BL_SE[af$site=="W1"], af$yr[af$site=="W1"]+0.15, af$BL_est[af$site=="W1"]+af$BL_SE[af$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 2002,
  y = min(af$BL_est-af$BL_SE-1) + 0.5,
  width = 2.5,
  height = NULL
)
title(main="b. Ancylus fluviatilis",bty="n",cex.main=2)

###############Baetis rhodani  
br[is.na(br)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(br$BL_est-br$BL_SE-1),max(br$BL_est + br$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$BL_est[br$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$BL_est[br$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$BL_est[br$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$BL_est[br$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$BL_est[br$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$BL_est[br$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$BL_est[br$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$BL_est[br$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(br$yr[br$site=="Auba"]-0.15, br$BL_est[br$site=="Auba"]-br$BL_SE[br$site=="Auba"], br$yr[br$site=="Auba"]-0.15, br$BL_est[br$site=="Auba"]+br$BL_SE[br$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="Bieb"]-0.05, br$BL_est[br$site=="Bieb"]-br$BL_SE[br$site=="Bieb"], br$yr[br$site=="Bieb"]-0.05, br$BL_est[br$site=="Bieb"]+br$BL_SE[br$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="O3"]+0.05, br$BL_est[br$site=="O3"]-br$BL_SE[br$site=="O3"], br$yr[br$site=="O3"]+0.05, br$BL_est[br$site=="O3"]+br$BL_SE[br$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="W1"]+0.15, br$BL_est[br$site=="W1"]-br$BL_SE[br$site=="W1"], br$yr[br$site=="W1"]+0.15, br$BL_est[br$site=="W1"]+br$BL_SE[br$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/baetis_rhodani.svg",
  x = 2002,
  y = min(br$BL_est-br$BL_SE-1) + 8.75,
  width = 4,
  height = NULL
)
title(main="c. Baetis rhodani",bty="n",cex.main=2)
            
###############Eiseniella tetraedra 
et[is.na(et)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(et$BL_est-et$BL_SE-1),max(et$BL_est + et$BL_SE+1)), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.7,cex.lab=2)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$BL_est[et$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$BL_est[et$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$BL_est[et$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$BL_est[et$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$BL_est[et$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$BL_est[et$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$BL_est[et$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$BL_est[et$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(et$yr[et$site=="Auba"]-0.15, et$BL_est[et$site=="Auba"]-et$BL_SE[et$site=="Auba"], et$yr[et$site=="Auba"]-0.15, et$BL_est[et$site=="Auba"]+et$BL_SE[et$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="Bieb"]-0.05, et$BL_est[et$site=="Bieb"]-et$BL_SE[et$site=="Bieb"], et$yr[et$site=="Bieb"]-0.05, et$BL_est[et$site=="Bieb"]+et$BL_SE[et$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="O3"]+0.05, et$BL_est[et$site=="O3"]-et$BL_SE[et$site=="O3"], et$yr[et$site=="O3"]+0.05, et$BL_est[et$site=="O3"]+et$BL_SE[et$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="W1"]+0.15, et$BL_est[et$site=="W1"]-et$BL_SE[et$site=="W1"], et$yr[et$site=="W1"]+0.15, et$BL_est[et$site=="W1"]+et$BL_SE[et$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/eiseniella_tetraedra.svg",
  x = 2002,
  y = min(et$BL_est-et$BL_SE-1) + 3.7,
  width = 3,
  height = NULL
)
title(main="d. Eiseniella tetraedra",bty="n",cex.main=2)
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)

###############Ephemera danica 
ed[is.na(ed)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ed$BL_est-ed$BL_SE-1),max(ed$BL_est + ed$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$BL_est[ed$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$BL_est[ed$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$BL_est[ed$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$BL_est[ed$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$BL_est[ed$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$BL_est[ed$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$BL_est[ed$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$BL_est[ed$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(ed$yr[ed$site=="Auba"]-0.15, ed$BL_est[ed$site=="Auba"]-ed$BL_SE[ed$site=="Auba"], ed$yr[ed$site=="Auba"]-0.15, ed$BL_est[ed$site=="Auba"]+ed$BL_SE[ed$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="Bieb"]-0.05, ed$BL_est[ed$site=="Bieb"]-ed$BL_SE[ed$site=="Bieb"], ed$yr[ed$site=="Bieb"]-0.05, ed$BL_est[ed$site=="Bieb"]+ed$BL_SE[ed$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="O3"]+0.05, ed$BL_est[ed$site=="O3"]-ed$BL_SE[ed$site=="O3"], ed$yr[ed$site=="O3"]+0.05, ed$BL_est[ed$site=="O3"]+ed$BL_SE[ed$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="W1"]+0.15, ed$BL_est[ed$site=="W1"]-ed$BL_SE[ed$site=="W1"], ed$yr[ed$site=="W1"]+0.15, ed$BL_est[ed$site=="W1"]+ed$BL_SE[ed$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/Ephemera_danica.svg",
  x = 2002,
  y = min(ed$BL_est-ed$BL_SE-1) + 2.75,
  width = 3.5,
  height = NULL
)
title(main="e. Ephemera danica ",bty="n",cex.main=2)

###############Gammarus roeselii
gr[is.na(gr)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(gr$BL_est-gr$BL_SE-1),max(gr$BL_est + gr$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$BL_est[gr$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$BL_est[gr$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$BL_est[gr$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$BL_est[gr$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$BL_est[gr$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$BL_est[gr$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$BL_est[gr$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$BL_est[gr$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(gr$yr[gr$site=="Auba"]-0.15, gr$BL_est[gr$site=="Auba"]-gr$BL_SE[gr$site=="Auba"], gr$yr[gr$site=="Auba"]-0.15, gr$BL_est[gr$site=="Auba"]+gr$BL_SE[gr$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="Bieb"]-0.05, gr$BL_est[gr$site=="Bieb"]-gr$BL_SE[gr$site=="Bieb"], gr$yr[gr$site=="Bieb"]-0.05, gr$BL_est[gr$site=="Bieb"]+gr$BL_SE[gr$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="O3"]+0.05, gr$BL_est[gr$site=="O3"]-gr$BL_SE[gr$site=="O3"], gr$yr[gr$site=="O3"]+0.05, gr$BL_est[gr$site=="O3"]+gr$BL_SE[gr$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="W1"]+0.15, gr$BL_est[gr$site=="W1"]-gr$BL_SE[gr$site=="W1"], gr$yr[gr$site=="W1"]+0.15, gr$BL_est[gr$site=="W1"]+gr$BL_SE[gr$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/Gammarus_roselii.svg",
  x = 2002,
  y = min(gr$BL_est-gr$BL_SE-1) + 1.85,
  width = 3,
  height = NULL
)
title(main="f. Gammarus roeselii",bty="n",cex.main=2)

###############Hydropsyche siltalai 
hs[is.na(hs)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(hs$BL_est-hs$BL_SE-1),max(hs$BL_est + hs$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$BL_est[hs$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$BL_est[hs$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$BL_est[hs$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$BL_est[hs$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$BL_est[hs$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$BL_est[hs$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$BL_est[hs$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$BL_est[hs$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(hs$yr[hs$site=="Auba"]-0.15, hs$BL_est[hs$site=="Auba"]-hs$BL_SE[hs$site=="Auba"], hs$yr[hs$site=="Auba"]-0.15, hs$BL_est[hs$site=="Auba"]+hs$BL_SE[hs$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="Bieb"]-0.05, hs$BL_est[hs$site=="Bieb"]-hs$BL_SE[hs$site=="Bieb"], hs$yr[hs$site=="Bieb"]-0.05, hs$BL_est[hs$site=="Bieb"]+hs$BL_SE[hs$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="O3"]+0.05, hs$BL_est[hs$site=="O3"]-hs$BL_SE[hs$site=="O3"], hs$yr[hs$site=="O3"]+0.05, hs$BL_est[hs$site=="O3"]+hs$BL_SE[hs$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="W1"]+0.15, hs$BL_est[hs$site=="W1"]-hs$BL_SE[hs$site=="W1"], hs$yr[hs$site=="W1"]+0.15, hs$BL_est[hs$site=="W1"]+hs$BL_SE[hs$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/hydropsyche_siltalai.svg",
  x = 2002,
  y = min(hs$BL_est-hs$BL_SE-1) + 1,
  width = 4,
  height = NULL
)
title(main="g. Hydropsyche siltalai",bty="n",cex.main=2)

###############Orectochilus villosus
ov[is.na(ov)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ov$BL_est-ov$BL_SE-1),max(ov$BL_est + ov$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$BL_est[ov$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$BL_est[ov$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$BL_est[ov$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$BL_est[ov$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$BL_est[ov$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$BL_est[ov$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$BL_est[ov$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$BL_est[ov$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(ov$yr[ov$site=="Auba"]-0.15, ov$BL_est[ov$site=="Auba"]-ov$BL_SE[ov$site=="Auba"], ov$yr[ov$site=="Auba"]-0.15, ov$BL_est[ov$site=="Auba"]+ov$BL_SE[ov$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="Bieb"]-0.05, ov$BL_est[ov$site=="Bieb"]-ov$BL_SE[ov$site=="Bieb"], ov$yr[ov$site=="Bieb"]-0.05, ov$BL_est[ov$site=="Bieb"]+ov$BL_SE[ov$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="O3"]+0.05, ov$BL_est[ov$site=="O3"]-ov$BL_SE[ov$site=="O3"], ov$yr[ov$site=="O3"]+0.05, ov$BL_est[ov$site=="O3"]+ov$BL_SE[ov$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="W1"]+0.15, ov$BL_est[ov$site=="W1"]-ov$BL_SE[ov$site=="W1"], ov$yr[ov$site=="W1"]+0.15, ov$BL_est[ov$site=="W1"]+ov$BL_SE[ov$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/Orectochilus_villosus.svg",
  x = 2002,
  y = min(ov$BL_est-ov$BL_SE-1) + 11.5,
  width = 2.5,
  height = NULL
)
title(main="h. Orectochilus villosus",bty="n",cex.main=2)

###############Prodiamesa olivacea  
po[is.na(po)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(po$BL_est-po$BL_SE-1),max(po$BL_est + po$BL_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$BL_est[po$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$BL_est[po$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$BL_est[po$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$BL_est[po$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$BL_est[po$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$BL_est[po$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$BL_est[po$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$BL_est[po$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(po$yr[po$site=="Auba"]-0.15, po$BL_est[po$site=="Auba"]-po$BL_SE[po$site=="Auba"], po$yr[po$site=="Auba"]-0.15, po$BL_est[po$site=="Auba"]+po$BL_SE[po$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="Bieb"]-0.05, po$BL_est[po$site=="Bieb"]-po$BL_SE[po$site=="Bieb"], po$yr[po$site=="Bieb"]-0.05, po$BL_est[po$site=="Bieb"]+po$BL_SE[po$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="O3"]+0.05, po$BL_est[po$site=="O3"]-po$BL_SE[po$site=="O3"], po$yr[po$site=="O3"]+0.05, po$BL_est[po$site=="O3"]+po$BL_SE[po$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="W1"]+0.15, po$BL_est[po$site=="W1"]-po$BL_SE[po$site=="W1"], po$yr[po$site=="W1"]+0.15, po$BL_est[po$site=="W1"]+po$BL_SE[po$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_phylopic(
  upload_img = "Silhouettes/Prodiamesa_olivacea.svg",
  x = 2002,
  y = min(po$BL_est-po$BL_SE-1) + 12,
  width = 3.5,
  height = NULL
)
title(main="i. Prodiamesa olivacea ",bty="n",cex.main=2)
dev.off()
##
##


