spann <- 0.2
require(stringr)
library(adespatial)
library(reshape)
library(plotrix)
library("vegan")
library(colorspace)
library(pgirmess)
library(dendextend)

library(analogue)
library(rioja)

colMax <- function (colData) { apply(colData, MARGIN=c(2), max)}
source('pollen.R')

### READING DATA, EXCLUDING NON TARGET TAXA 
    si <- read.table( "cl4_pollen_body_mapa.txt")
    ag <- read.table("agde.txt")
    n <- read.table("numbers.txt", stringsAsFactors = F)
    
    
    n <- n[!(n$anthropo %in% c("Cyperaceae", "aaax", "NPP")),]
    
    lok <- c("Butomus umbellatus"  ,   "Alnus",     "Hippuris vulgaris"          ,
     "Hydrocotyle vulgaris"  ,      "Alisma-Typ"                 ,
     "Litorella-Typ"          ,     "Lobelia dortmanna"          ,
     "Menyanthes trifoliata"   ,    "Myriophyllum alterniflorum" ,
     "Myriophyllum spicatum"    ,   "Myriophyllum verticillatum" ,
     "Nuphar"                    ,  "Nymphaea"                   ,
     "Varia"                      , "Lemnaceae"                  ,
     "Potamogeton natans-Typ"      ,"Sagittaria sagittifolia"    ,
     "Scheuchzeria palustris",      "Potamogeton pectinatus-Typ" ,
     "Sparganium-Typ"         ,     "Stratiotes aloides"         ,
     "Typha latifolia-Typ"      ,   "Ludwigia palustris"         ,
     "Cladium"                 ,    "Calla palustris"            ,
     "Caltha-Typ",                  "Cyperaceae"                 ,
     "Drosera"              ,
     "Iris pseudacorus-Typ" ,       "Berula erecta, Sium sisarum",
     "Utricularia"            ,     "Callitriche"                ,
     "Trapa natans"          ,      "Triglochin"                 ,
     "Hottonia palustris"      ,    "Aldrovanda vesiculosa"      ,
     "Hydrocharis morsus-ranae" ,   "Nymphoides peltata"       ,
     
     "Ranunculaceae indet",    "Primulaceae",  "Plantaginaceae", "Tofieldia" 
     )
    
    n <- n[!(n$anthropo %in% lok),]
    
    n[n$anthropo == "Sammelgruppe Apiaceae","anthropo"] <- "Apiaceae"
    n[n$anthropo == "Sammelgruppe Vaccinium-Typ","anthropo"] <- "Vaccinium-Typ"   
    
    
    trees <- c("Sorbus/Cotoneaster" ,    "Taxus"  ,             "Picea",                 
     "Carpinus betulus"   ,    "Ulmus"                 , "Fraxinus excelsior-Typ",
     "Acer"                ,   "Tilia"              ,    "Quercus"   ,            
     "Pinus"              ,    "Alnus"               ,   "Betula"   ,             
                   "Fraxinus ornus"   ,     
     "Abies"                ,  "Corylus",                "Populus",               
     "Fagus")               
    
    
    druhy <- c("Pinus" , "Betula"  ,"Corylus",  "Picea"  , 
               "Quercus" ,   "Ulmus",      "Fraxinus excelsior-Typ" ,  "Tilia",                 
               "Abies", "Fagus",               "Carpinus betulus","Wildgras-Typ",  
               "Cerealia", "Secale","Plantago lanceolata-Typ", "Chenopodiaceae pp", "Artemisia","Compositae Subfam Cichorioideae",
               "Polygonum aviculare-Typ" ,  "Centaurea cyanus type", "Urtica"
               
    )
    #aunique(n$anthropo)
    
    n <- n[(n$anthropo %in% druhy),]
    nt  <- aggregate(n$count, by=list(n$e, n$depth, n$anthropo), FUN=sum)  # get rid of precision of individual pollen analysts 
    colnames(nt) <- colnames(n)[-3]

### MAKING TIME WINDOWS
        sm <- aggregate(nt$count, by=list(nt$e, nt$depth), FUN=sum) 
        sma <- merge(ag[,c(1,3,4)], sm, by=c(1,2), all.y=T)
        sma <- sma[!is.na(sma$age),]
        sma <- sma[order(sma$e., sma$depth),]
        di <- data.frame(dyf=NA, site=NA)
        di <- di[-1,]
        i=1
        for (i in 1:(NROW(si))) {
          dyf <- diff(sma[sma$e.==si[i,1], "age"])
          di <- rbind(di, data.frame(dyf=dyf, site=rep(si[i,"Idword"], NROW(dyf))) )
                                                                 }
        di[is.na(di$dyf),]
        
        
        delka_obdobi <- 500
        obdobi <- seq(-11000,2000,delka_obdobi) 
        obdobi <-  -1*(obdobi-1950)
        mids <- obdobi[1:(length(obdobi)-1)]+delka_obdobi/2
        cas <- cbind(1:(length(obdobi)-1), rev(obdobi)[-1], rev(obdobi)[-length(obdobi)]  )
        cs <- data.frame(cas[,c(1,3,2)])
        
        i=1
        sma[,5] <- paste(sma$e., sma$depth, sep="_")
        smtw <- age_to_tw(sma, "age", cs, 2,3)
        smtw[,7] <- paste(smtw$e., smtw$tw, sep="_") 
        smtw[,8] <- smtw$x
        smtwa <- smtw
        
        nt[, ncol(nt)+1]<- paste(nt$e, nt$depth, sep="_")
       
### AVERAGING WITHING TW & REGIONS  
        nrtw <- merge(nt, smtwa[,c( "V5", "tw", "x")], by.x=5,by.y=1, all.y = T)
        nrtw <- aggregate(nrtw[,c("count", "x")], by=list(nrtw$e, nrtw$tw, nrtw$anthropo), FUN=sum)
        nrcl <- merge(nrtw, si[,c("ID", "class10")], by.x=1, by.y=1)
        obili <- c("Zea mays", "Secale", "Cerealia")
    
        nrca <- aggregate(nrcl[,c("count","x")], by=list( nrcl$class10, nrcl$Group.2,nrcl$Group.3), FUN=sum)
        nrca$pp <- nrca$count/nrca$x
        pct <- unique(nrcl[,c(1:2,6)])
        pct <- aggregate(pct$Group.1, by=list(pct$class10, pct$Group.2), FUN=length)
        pct <- cast(pct, Group.2~Group.1)
  
### READING ARCHAEOLOGY
        archmean <- read.csv("class_mean.csv", sep=";", check.names = F)
        archm <- data.frame(age=(seq(-9250, 1250, 500)-1950)*-1,t(archmean[,c(6:27)]))
        #archsd <- read.csv("class_sd.csv", sep=";", check.names = F)
        #archs <- data.frame(age=(seq(-9750, 1750, 500)-1950)*-1,t(archsd[,c(5:28)]))
        colnames(nrca)[1:3] <- c("reg", "age", "variable")
        nn <- unique(nrca$reg)
        tt <- unique(nrca$age)
        tt <- tt[tt!=200]
      
### PLOTING DATA OVERVIEW        
        names<- c("sparsely occ. area", "Plzen", "W Polabi", "South Bohemia", "E Polabi", "S Moravia", "Boskovice" , "Hornomoravsky", "White Carp.", "Opava")
        druhy <- c("Pinus" , "Betula"  ,"Corylus",  "Picea"  , 
        "Quercus" ,   "Ulmus",      "Fraxinus excelsior-Typ" ,  "Tilia",                 
        "Abies", "Fagus",               "Carpinus betulus","Wildgras-Typ",  
         "Cerealia", "Secale","Plantago lanceolata-Typ", "Chenopodiaceae pp"
        )
        drh <- c("Pinus" , "Betula"  ,"Corylus",  "Picea"  , 
                   "Quercus" ,   "Ulmus",      "Fraxinus" ,  "Tilia",                 
                   "Abies", "Fagus",               "Carpinus","wild grasses",  
                   "Cerealia", "Secale","Plantago lanceolata", "Chenopodiaceae", "occupancy"
        )
        
        fnt <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4 ,4 ,2 ,2 ,4 ,4 ,2, 2) 
        colc <- c("#6AB2A6", "#f3c300" ,"#2b3d26", "#DD6555", "#6295B7", "#DA9437", "#91BB3C", "#D7A9C1", "gray40", "#A467A5")
        names<- c("Sparse o.a.", "Plzen", "C & NW B.", "S Bohem.", "E Polabi", "S Morav.", "Boskovice" , "C Morav.", "White Carp.", "Opava")
        vyb <- cast(nrca[nrca$reg %in% c(1,3,4,6,8),],formula=reg+age~variable) 
        vyb[is.na(vyb)] <- 0
        vyb <- vyb[,c("reg", "age",druhy) ]
        
        aroc <- melt(archm,id="age")
        aroc <- data.frame(Group.1=as.integer(gsub("X","", as.character(aroc$variable))),
                           aroc[,-2])
        vyb <- merge(vyb, aroc, by=1:2)#, all.x = T)
        tt <- rev(seq(700,11200,500))
        popi <- (abs((tt+250)-1950)[seq(2,NROW(tt), 2)])
        popi[10] <- 1
        
  tiff("pollen_averaged.tif",3500,1800,res=300,comp="lzw", pointsize = 12)
        par(mar=c(7,0.5,20,0), mfrow=c(1,19), oma=c(0,2,0,2))
        y <- rev(sort(tt))
        
        plot(1:2,1:2,  col= 'transparent',  ylab="", xlab="", xlim=c(0, 20 ), ylim=c(y[NROW(y)],y[1]), xaxt="n",yaxt="n",xaxt="n", main="", bty="n")
        for (i in 1:(NROW(druhy)+1))#(length(jmn)))
        {
          mxx <- max(vyb[, i+2])*100
          plot(rep(0,NROW(y)), y, "l",ylim=c(y[1], y[NROW(y)]), xlab="", ylab="", xlim=c(0,mxx), yaxt="n",xaxt="s",lwd=3, bty="n", col="transparent")
          abline(h=c(8200,4200,350), col="gray", xpd=T, lty=2, lwd=2)
          text(cex=2, x=0, y=-500+5, drh[i], xpd=TRUE, srt=90, pos=4, font=fnt[i] )
         v=1
          for(v in c(1,3,4,6,8)){
            vybv <- merge(tt, vyb[vyb$reg==v,-1], by=1, all.x=T)
            vybv <- vybv[order(-vybv$x),-1 ]*100#c(1,3,2,4)]
            lines(vybv[,i],y, col=colc[v], lwd=2)  
           }
          
          if(i==1){          axis(2, at=seq(950,10950, 1000), labels=rev(popi), las=1, cex.axis=1.5)      }   
          if(i==(NROW(druhy)+1)){
            axis(4, at=c(9850, 6200, 2200, 100), labels=c("Early H.", "Middle H.", "Late H.", "Modern P."), las=1, cex.axis=1.5, tick = F)
            legend("bottom", ncol = 1,pch=15, col=colc[c(1,3,4,6,8)], 
                   legend = names[c(1,3,4,6,8)], xpd = TRUE,  inset = c(0, -0.38), bty = "n", cex = 0.9)
            
            }  
        }
        mtext("age (BCE/CE)",outer=T,side=2,cex=1.3, font=1, line=0, adj=0.3)
        mtext("%" ,outer=T,side=1,cex=1.3, font=1, line=-2, adj=0.5)
        dev.off()
        
        
### CLUSTER & BETA DIVERSITY ANALYSIS
      age <- tt
      druhy <- c("Pinus" , "Betula"  ,"Corylus",  "Picea"  , 
                 "Quercus" ,   "Ulmus",      "Fraxinus excelsior-Typ" ,  "Tilia",                 
                 "Abies", "Fagus",               "Carpinus betulus","Wildgras-Typ",  
                 "Cerealia", "Secale","Plantago lanceolata-Typ", "Chenopodiaceae pp", "Artemisia","Compositae Subfam Cichorioideae",
                 "Polygonum aviculare-Typ" ,  "Centaurea cyanus type", "Urtica"
                 
      )
      
        nnam <- names[c( 3,4,6, 8, 1)]
        vyb <- nrca[nrca$variable %in% druhy & nrca$reg %in% c(1,3,4,6, 8) & nrca$age %in% age,]#c(1,2,3,6)]
        vyb <- cast(vyb, formula = reg+age~variable)
        vyb[is.na(vyb)] <-0
        prumz <- vyb
        rownames(prumz) <- paste(prumz[,1], prumz[,2], sep="_")
        prumzonREV <- prumz[,-1]
        prumzonREV <- prumzonREV[,-1]
        
        prumzonREV <- pclig(prumzonREV)
        sqd <- hclust(as.dist(as.matrix(distance(prumzonREV, method="SQchord")[1:107, 1:107])), "ward.D")
        clust <- cutree(sqd, 6)
        clust.cutree <- dendextend:::cutree(sqd, k=6, order_clusters_as_data = FALSE)
        idx <- order(names(clust.cutree))
        clust.cutree <- clust.cutree[idx]
        tbl <- table(clust, clust.cutree)
        lbls <- apply(tbl,2,which.max)
        lbls[6]<-2
        
        naz <- merge( as.data.frame(cutree(sqd,6)),prumzonREV, by=0)[,1:2]
        prumz[,24] <- (paste(prumz[,1], prumz[,2], sep="_"))
        nbx <- merge(prumz, naz, by.x=24, by.y=1)
        
        doimg <- merge(prumz[,c(24,1,2)], naz, by=1, sort=F, all=T)
        
        
        
        bar<- c("#ff3700", # cervena
                "#ff9509", # oranzova vyblita
                "#084086" , # tm modra
                "lightgreen", # sv zelena 
                "#89358a", # magenta
                "#ffe12b" )#,#zluta vyblita
        
        imn <- cast(doimg,reg~age)
        rownames(imn) <- imn[,1]
        imn <- imn[,-1]
        imn <- imn[ c("3","4","6", "8", "1"),rev(colnames(imn))]
        im <- imn
        
        
        nadendro <- list(imn, sqd)
        imn <- nadendro[[1]]
        im <- imn
        sqd <- nadendro[[2]]
        
        
        
        vari_spat <- numeric()
        vari_hol <- numeric()
        vari_arch <- numeric()
        arch_spat <- archm[order(archm$age), c("X3","X4","X6", "X8")]
        arch_arch <- archm[order(archm$age), ]
        arch_arch <- arch_arch[arch_arch$age<=7200, c("X3","X4","X6", "X8")]
        for (i in 1:nrow(arch_spat)) {vari_spat[i] <- var(t(arch_spat)[,i])}
        for (i in 1:ncol(arch_spat))  {vari_hol[i]  <- var(t(arch_spat)[i,])}
        for (i in 1:ncol(arch_arch))  {vari_arch[i]  <- var(t(arch_arch)[i,])}
        
        
        
        vybw <- vyb[vyb$reg!=1,]
        bet_p <- numeric()
        drh <- colnames(vybw)[3:ncol(vybw)]
        sc_p <- data.frame(matrix(nrow=NROW(tt),ncol=NROW(drh)), row.names = tt) 
        colnames(sc_p) <- drh
        lc_p  <- data.frame(matrix(nrow=0,ncol=4)) 
        colnames(lc_p) <- c("regid","tw" ,"lcp", "plcp")
        me = "chord"
            
        for(i in NROW(tt):1){
          res <- beta.div(vybw[vybw$age==tt[i],3:ncol(vybw)], method = me)
          bet_p[i] <- res$beta["BDtotal"]
          sc_p[i,] <-  res$SCBD
          lc_p <-  rbind(lc_p,data.frame(regid=vybw[vybw$age==tt[i],1] , tw=rep(tt[i],NROW(res$LCBD)),res$LCBD,res$p.adj)) 
        }  
        betres <- list(bet_p,sc_p,lc_p, vari_spat)
        
        bet_p <- betres[[1]]
        vari <- betres[[4]]
        cor.test(rev(bet_p), vari)
        cor.test(rev(bet_p[1:14]), vari[1:14])
        lcp <- betres[[3]][,1:3]
        lcp[]
        require(reshape)
        lcp <- cast(lcp, formula=tw~regid)
        rownames(lcp) <- lcp[,1]
        lcp <- lcp[order(-lcp$tw),-1]
 
### ANALYSIS FINISHED
### PLOTTING:
        
        
  tiff("dendro2.tiff",width = 400, height = 300, units="mm", res=300, pointsize=25, compression="lzw" )
        par(mar=c(1,5,2,1), mfrow=c(1,2))
        sqd %>%
          color_branches(sqd,k=6, col = bar[  c(3,5, 4,6,1,2)   ],  groupLabels = lbls) %>%# auto-coloring 5 clusters of branches.
          set("branches_lwd", c(6))   %>%                  #kodskup[  3,4,2,1, 5,6   ]
          plot(cex=0.1, xlab=" ", ylab="SQchord distance", main="", las=1, leaflab = "none")
        mtext("a)", side = 3, outer = F, line=0, adj=0, cex = 1.5)
        par(mar=c(2,5,10,4))
        image(as.matrix(im), col=bar, axes=F, ylab="age (BC/AD)")
        mtext("b)", side = 3, outer = F, line=8, adj=0, cex = 1.5)
        axis(2, at=rev(seq(0,1,1/(nrow(lcp)-1))),abs(as.integer(rev(rownames(lcp)))-1950) , las=2)
        axis(3, at=    seq(0,1,1/(nrow(im)-1)), nnam, cex.axis=1.2, line=0, tick=F, las=2)
        #axis(2, at=rev(seq(0,1,1/(ncol(im)-1)) -((1/(ncol(im)-1))/2))[seq(1,ncol(im),4)],seq(400,15400,500)[seq(1,ncol(im),4)] , las=2, cex.axis=.9)
        #abline(h=  rev(seq(0,1,1/(ncol(im)-1)) -((1/(ncol(im)-1))/2))[seq(1,ncol(im),4)], col="black", xpd=F, lty=2, lwd=1)
        dev.off()
        
        
  tiff("boxplot.tiff",width = 650, height = 350, units="mm", res=300, pointsize=22, compression="lzw" )
        par(mfrow=c(4,6))
        for(i in 3+1:21){
          boxplot(nbx[,i]~nbx[,25], main=colnames(nbx)[i], ylab="prop.", xlab= "class", col=bar)
        }
        dev.off()
        

        lbls <- c(3,5,4,6,1,2)
        sloupcu=5
        sir=c(0.6,0.5,0.5,0.8,0.9,1.8,0.3)
        tt <- as.integer(colnames(im))
 tiff("envi_analw.tiff",width = 500, height = 300, units="mm", res=300, pointsize=37, compression="lzw" )
        nf <- layout(matrix( seq(1,1*sloupcu+2),1,sloupcu+2,byrow=T), width=sir, height=c(5))
        layout.show(nf)
        par(mar=c(4,0,12,1) , pty="m")
        y <- barplot(rev(vari), space = 0, border=F,col= 'transparent', horiz = T,  ylab="", xlab="",  xaxt="n",yaxt="n",xaxt="n", main="", bty="n",  xaxs="i",yaxs="i")
        mtext("a)", side = 3, outer = F, line=10, adj=0.5, cex = 1.1)
        
        barplot(rev(vari), horiz=T, las=2, space = 0, border=F,col= 'black',  ylab="", xlab="",  yaxt="n", main="", bty="n",  xaxs="i",yaxs="i")
        popi <- (abs((tt+250)-1950)[seq(2,NROW(tt), 2)])
        popi[10] <- 1
        axis(2, at=y[seq(2,NROW(y), 2)]-0.5,labels = popi  , las=2,  xpd=T,  mgp=c(1.5,0.5,0))
        axis(3,line=0, at=max(vari)/6,labels="variance of occupancy",tick=F, las=2)
        
        mtext("age (BCE/CE)", side = 2, outer = T, line=-1.5, adj=0.4, cex = 0.8)
        
        barplot((bet_p), horiz=T,  main="", las=2, space = 0, border=F,col= 'black',xlab = "", xpd=F, las=2,  xaxs="i",yaxs="i")
        axis(3, line=0,at=max(bet_p)/6,labels="variance of pollen",tick=F, las=2)
        mtext("b)", side = 3, outer = F, line=10, adj=0, cex = 1.1)
        
        image(t(as.matrix(lcp)),col=rev(heat.colors(12)), axes=F, ylab="age cal. BP", xlab="")
        betce <- seq(0,1,1/(ncol(lcp)-1))
        axis(3, at= c(  betce ),names[as.integer(colnames(lcp))],  line=0, tick=F, las=2)
        mtext("c)", side = 3, outer = F, line=10, adj=0, cex = 1.1)
        
        legend("bottom",  fill=c("yellow", "red"), legend = c("low", "high"), xpd = TRUE, horiz = F, inset = c(0, -0.25), bty = "n")
        
        
        image(as.matrix(im), col=bar, axes=F, ylab="Age (years BP)")
        axis(3, at=    seq(0,1,1/(nrow(im)-1)), nnam,  line=0, tick=F, las=2)
        abline(v=0.875, lty=2, lwd=3)
        axis(4, at=rev(seq(0,1,1/(nrow(lcp)-1))+(1/(nrow(lcp)-1))/2 )[seq(2,nrow(lcp), 2)],rev(popi) , las=2,  mgp=c(1.5,0.5,0))
        mtext("d)", side = 3, outer = F, line=10, adj=0, cex = 1.1)
        #axis(2, at=rev(seq(0,1,1/(ncol(im)-1)) -((1/(ncol(im)-1))/2))[seq(1,ncol(im),4)],seq(400,15400,500)[seq(1,ncol(im),4)] , las=2, cex.axis=.9)
        #abline(h=  rev(seq(0,1,1/(ncol(im)-1)) -((1/(ncol(im)-1))/2))[seq(1,ncol(im),4)], col="black", xpd=F, lty=2, lwd=1)
        legend("bottom", ncol = 3,pch=15, col=bar, legend = 1:6, xpd = TRUE,  inset = c(0, -0.25), bty = "n")
        
        
        par(mar=c(3,6,3,1), mgp=c(1.5,0.5,0))
        sqd %>%
          color_branches(sqd,k=6, col = bar[  c(3,5, 4,6,1,2)   ],  groupLabels = lbls) %>%# auto-coloring 5 clusters of branches.
          set("branches_lwd", c(6))   %>%                  #kodskup[  3,4,2,1, 5,6   ]
          plot(cex=0.1, xlab=" ", ylab="SQchord distance", main="", las=1, leaflab = "none")
        mtext("e)", side = 3, outer = F, line=1, adj=0, cex = 1.1)
        
        legend("bottom", ncol = 3,pch=15, col=bar, legend = 1:6, xpd = TRUE,  inset = c(0, -0.1), bty = "n")
        dev.off()
