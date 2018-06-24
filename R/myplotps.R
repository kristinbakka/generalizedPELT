## Functions for plot and save
## From 21c03.R
# Function for plot
myplot <- function(data,name,runs,ylab="Total cost per length",newframe=FALSE){
  #myplot(data=l.plot[[i]],name=names(l.plot)[i],runs)
  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$variable))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }
  return(ggplot(data=data,
                #aes(x=n, y=value, group=variable,colour=variable)) +
                aes(x=n, y=value, colour=variable)) +
           geom_point(aes(group = variable))+
           geom_line(aes(group = variable))+theme_minimal()+
           #coord_cartesian(ylim=c(0,1.10))+
           #scale_fill_continuous()+
           scale_color_manual("Method",values=mypal)+
           guides(color=guide_legend(title="penalty"))+
           labs(x="n",y=ylab,title=paste("runs",runs,name))
  )
}

barplot.single <- function(data,name,runs,ylab="Total cost per length",newframe=FALSE){
  data$delta=as.character(data$delta)


  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$delta))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data, aes(x=ncpts,y=Freq,
                                          color=delta,group=delta)) +
    #geom_point(position=position_dodge(width = 0.9))+
    geom_point()+
    geom_line(linetype="dotted")+
    theme_minimal()+
    scale_color_manual(values=mypal)+
    #  scale_colour_gradient()+
    labs(x="m",colour=expression(Delta),title=paste("runs",runs,name),y = "Proportion", fill = "Delta")+
    #coord_flip(ylim = c(0,1))
    coord_cartesian(ylim=c(0,1))
  )
}

barplot2.single <-function(data,name,runs,ylab="Total cost per length",newframe=FALSE){
  data$ncpts=as.character(data$ncpts)
  data$delta=as.numeric(data$delta)

  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$ncpts))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data, aes(x=delta,y=Freq,
                               color=ncpts,group=ncpts)) +
           #geom_point(position=position_dodge(width = 0.9))+
           geom_point()+
           geom_line()+
           theme_minimal()+
           scale_color_manual(values=mypal)+
           #  scale_colour_gradient()+
           labs(x=expression(Delta),colour="m",title=paste("runs",runs,name),y = "Proportion", fill = "Delta")+
           #coord_flip(ylim = c(0,1))
           coord_cartesian(ylim=c(0,1))
  )
}

barplot2.double <-function(data,name,runs,newframe=FALSE){
  data$ncpts=as.character(data$ncpts)
  data$delta=as.numeric(data$delta)

  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$ncpts))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data, aes(x=delta,y=Freq,
                               color=ncpts,group=interaction(ncpts,Method),
                               shape=Method, linetype=Method)) +
           #geom_point(position=position_dodge(width = 0.9))+
           scale_shape_manual(name="Method",values=c(4,19))+
           geom_point()+
           geom_line()+
           theme_minimal()+
           scale_color_manual(values=mypal)+
           #scale_color_discrete()+
           #  scale_linetype_discrete()+
           scale_linetype_manual(values=c(3,1))+
           #  scale_x_continuous(breaks = as.integer(dataset.breaks))+
           #  scale_x_discrete(breaks = dataset.breaks)+
           #scale_shape()+
           #  scale_colour_gradient()+
           theme(legend.spacing.y = unit(0, "cm"))+
           theme(legend.margin=margin(t = 0, unit='cm'))+
           labs(x=expression(Delta),colour="m",title=paste("runs",runs,name),y = "Proportion", fill = "Delta")+
           #coord_flip(ylim = c(0,1))
           coord_cartesian(ylim=c(0,1))

  )
}


barplot2.several <-function(data,name,runs,newframe=FALSE){
  data$ncpts=as.character(data$ncpts)
  data$delta=as.numeric(data$delta)

  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$ncpts))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data, aes(x=delta,y=Freq,
                               color=ncpts,group=interaction(ncpts,Method),
                               shape=Method, linetype=Method)) +
           #geom_point(position=position_dodge(width = 0.9))+
           scale_shape(name="Method")+
           geom_point()+
           geom_line()+
           theme_minimal()+
           scale_color_manual(values=mypal)+
           #scale_color_discrete()+
           #  scale_linetype_discrete()+
           scale_linetype()+
           #  scale_x_continuous(breaks = as.integer(dataset.breaks))+
           #  scale_x_discrete(breaks = dataset.breaks)+
           #scale_shape()+
           #  scale_colour_gradient()+
           theme(legend.spacing.y = unit(0, "cm"))+
           theme(legend.margin=margin(t = 0, unit='cm'))+
           labs(x=expression(Delta),colour="m",title=paste("runs",runs,name),y = "Proportion", fill = "Delta")+
           #coord_flip(ylim = c(0,1))
           coord_cartesian(ylim=c(0,1))

  )
}

barplot3.several <-function(data,name,runs,newframe=FALSE){
  data$ncpts=as.character(data$ncpts)
  data$delta=as.numeric(data$delta)

  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$Method))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data, aes(x=delta,y=Freq,
                               color=Method,group=interaction(ncpts,Method),
                               shape=Method, linetype=ncpts)) +
           #geom_point(position=position_dodge(width = 0.9))+
           scale_shape(name="Method")+
           geom_point()+
           geom_line()+
           theme_minimal()+
           scale_color_manual(name="Method",values=mypal)+
           #scale_color_discrete()+
           #  scale_linetype_discrete()+
           scale_linetype(name="m")+
           #  scale_x_continuous(breaks = as.integer(dataset.breaks))+
           #  scale_x_discrete(breaks = dataset.breaks)+
           #scale_shape()+
           #  scale_colour_gradient()+
           theme(legend.spacing.y = unit(0, "cm"))+
           theme(legend.margin=margin(t = 0, unit='cm'))+
           labs(x=expression(Delta),colour="m",title=paste("runs",runs,name),y = "Proportion", fill = "Delta")+
           #coord_flip(ylim = c(0,1))
           coord_cartesian(ylim=c(0,1))

  )
}


hitprop.plot.v01 <-function(data,name,runs,newframe=FALSE,stip=FALSE){
  if(newframe){set.window(FALSE)}
  n.line=2

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  if(stip){
    return(ggplot(data=data,
                aes(x=delta, y=Value, colour=Type,shape=Method,linetype=Type)) +
           geom_line()+theme_minimal()+
           geom_point()+
           scale_color_manual(
             name="Type",
             values=wes_palette( name="Cavalcanti"),
             labels=c("Exact", "Correct m"))+
           scale_linetype_discrete(
             name="Type",
             labels=c("Exact", "Correct m"))+
           scale_shape_manual(name="Method",
                              values=c(19,4))+
           ylim(0,1)+
           labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
    )}
  return(ggplot(data=data,
                aes(x=delta, y=Value, colour=Type,shape=Method)) +
           geom_line()+theme_minimal()+
           geom_point()+
           scale_color_manual(
             name="Type",
             values=wes_palette( name="Cavalcanti"),
             labels=c("Exact", "Correct m"))+
           scale_shape_manual(name="Method",
                              values=c(19,4))+
           ylim(0,1)+
           labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
  )
}


hitprop.plot.v01.multiple <-function(data,name,runs,newframe=FALSE,stip=FALSE){
  if(newframe){set.window(FALSE)}
  n.line=2

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  return(ggplot(data=data,
                aes(x=delta, y=Value, colour=Type,shape=Method)) +
           geom_line()+theme_minimal()+
           geom_point()+
           scale_color_manual(
             name="Type",
             values=wes_palette( name="Cavalcanti"),
             labels=c("Exact", "Correct m"))+
           scale_shape_discrete(name="Method",
                              solid=TRUE)+
           ylim(0,1)+
           labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
  )
}


hitprop.plot.v02 <-function(data,name,runs,newframe=FALSE,stip=FALSE){
  if(newframe){set.window(FALSE)}
  n.line=2

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  if(stip){
    return(ggplot(data=data,
                  aes(x=delta, y=Value, colour=Type,shape=Method,linetype=Type)) +
             geom_line()+theme_minimal()+
             geom_point()+
             scale_color_manual(
               name="Type",
               values=wes_palette( name="Cavalcanti"),
               labels=c("Correct m"))+
             scale_linetype_discrete(
               name="Type",
               labels=c("Correct m"))+
             scale_shape_manual(name="Method",
                                values=c(19,4))+
             ylim(0,1)+
             labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
    )}
  return(ggplot(data=data,
                aes(x=delta, y=Value, colour=Type,shape=Method)) +
           geom_line()+theme_minimal()+
           geom_point()+
           scale_color_manual(
             name="Type",
             values=wes_palette( name="Cavalcanti"),
             labels=c("Correct m"))+
           scale_shape_manual(name="Method",
                              values=c(19,4))+
           ylim(0,1)+
           labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
  )
}


hitprop.plot <-function(data,name,runs,newframe=FALSE,stip=FALSE){
  ## This is wrong
  warning("This is the wrong plot, I got confused when I tried to make haste <3.")
  if(newframe){set.window(FALSE)}
  n.line=2

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }

  if(stip){
    return(ggplot(data=data,
                  aes(x=delta, y=Value, colour=Type,shape=Method,linetype=Type)) +
             geom_line()+theme_minimal()+
             geom_point()+
             scale_color_manual(
               name="Method",
               values=wes_palette( name="Cavalcanti"),
               labels=c("PELT","BinSeg"))+
             scale_linetype_discrete(
               name="Method",
               values=wes_palette( name="Cavalcanti"),
               labels=c("PELT","BinSeg"))+
             scale_shape_manual(name="Type",
                                values=c(19,4),
                                labels=c("Exact", "Correct m"))+
             ylim(0,1)+
             labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
    )}
  return(ggplot(data=data,
                aes(x=delta, y=Value, colour=Type,shape=Method)) +
           geom_line()+theme_minimal()+
           geom_point()+
           scale_color_manual(
             name="Method",
             values=wes_palette( name="Cavalcanti"),
             labels=c("PELT","BinSeg"))+
           scale_shape_manual(name="Type",
                              values=c(19,4),
                              labels=c("Exact", "Correct m"))+
           ylim(0,1)+
           labs(x=expression(Delta),y="Proportion",title=paste("runs",runs,name))
  )
}


myplot.m <- function(data,name="",runs,ylab="Hit prop",newframe=FALSE){
  data$n=as.character(data$n)
  data$n=as.integer(data$n)
  #myplot(data=l.plot[[i]],name=names(l.plot)[i],runs)
  if(newframe){set.window(FALSE)}
  n.line=length(unique(data$m))

  if(n.line>4){
    mypal=wes_palette("Cavalcanti", n.line, type = "continuous")
  }else{
    mypal=wes_palette( name="Cavalcanti")
  }
  return(
    ggplot(data=data,
                #aes(x=n, y=value, group=n,colour=n)) +
                aes(x=m, y=value, colour=n)) +theme_minimal()+
           geom_point(aes(group = n))+
           #geom_line(aes(group = n))+
           #coord_cartesian(ylim=c(0,1.10))+
           #scale_color_manual("Method",values=mypal)+
           scale_fill_continuous()+
           #guides(color=guide_legend(title="n"))+
           labs(x="m",y=ylab,title=paste("runs",runs,name))
  )
}

# Function for plot
myplot.vars.g <- function(data,name,runs){
  #myplot(data=l.plot[[i]],name=names(l.plot)[i],runs)
  data$variable=as.integer(data$variable)
  data$g=as.character(data$g)
  #colnames(data)[2]="pen"

  return(ggplot(data=data,
                aes(x=variable, y=value, group=g, colour=g)) +
           geom_point()+geom_line()+theme_minimal()+
           #coord_cartesian(ylim=c(0,1.10))+
           #             scale_color_manual("Method",values=PuBu)+
           guides(color=guide_legend(title="g"))+
           labs(x="penalty",y="Variability in total cost per length",title=paste("runs",runs,name))
  )
}

# Function for plot
myplot.vars <- function(data,name,runs){
  #myplot(data=l.plot[[i]],name=names(l.plot)[i],runs)
  return(ggplot(data=data,
                aes(x=pen, y=value)) +
           geom_point()+geom_line()+theme_minimal()+
           #coord_cartesian(ylim=c(0,1.10))+
           #  scale_color_manual("Method",values=wes_palette( name="Cavalcanti"))+
           #guides(color=guide_legend(title="nnodes"))+
           labs(x="pen",y="Variability in total cost per length",title=paste("runs",runs,name))
  )
}

# Function for save
mysave <- function(plot,variety, save=FALSE,n,runs,filenam,m=0,wd=0){
  setwd("C:/Users/Kristin/Documents/v18-master/pelt-package-v01/plot")
  name=paste0(filenam,"-m",m,"-n",n,"-r",runs,"-",variety) #i.e. "p21-r50-mean","p21-r50-meanvar"

  if(save){
    ggsave(plot=plot,  paste(c(name,".pdf"), collapse = ""), device="pdf")
    return("The pdf was saved.")
  }
  return("The pdf was not saved.")
}

mysave.proj <- function(plot,filenam,m=0,runs,variety,folder="plot-06-part2/session0",save=FALSE,wd=TRUE){
  # opens new window, plots in it, shuts it
       #-(what is now the active window?)
  set.window(FALSE)
  my.plot<- plot + labs(subtitle=paste0("m= ",m,", runs=",runs,", variety=",variety))
  clean.plot<- plot + ggtitle("")
  workdir="C:/Users/Kristin/Documents/v18-master/pelt-package-v01/"
  if(wd){setwd(paste0(workdir,folder))}
  name=paste0(filenam,"-m",m,"-r",runs,"-",variety) #i.e. "p21-r50-mean","p21-r50-meanvar"

  if(save){
    ggsave(plot=plot,  paste0(name,".pdf"), device="pdf", width=5)
    ggsave(plot=plot+ggtitle(""),  paste0("clean-",name,".pdf"), width=5, device="pdf")
    dev.off()
    return("The pdf was saved.")
  }
  return("The pdf was not saved.")
}



