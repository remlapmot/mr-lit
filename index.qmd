library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(ggplot2)

## plot settings
theme_mrlit<-function(...){
  theme_classic(base_size=14)%+replace%
    theme(legend.position="bottom",
          strip.background=element_rect(colour='grey99',fill='grey98'),
          panel.grid.major=element_line(linewidth=lwidth/3,colour='grey'),
          panel.grid.minor=element_line(linewidth=lwidth/5,colour='grey'),
          plot.background=element_rect(fill=NA,linewidth=lwidth/2,colour='grey'),
          panel.border=element_rect(fill=NA,linewidth=lwidth/2,colour='grey'),
          axis.line=element_line(linewidth=.2),
          axis.ticks=element_line(linewidth=.2))
}

## static vars
fpath="."
fname="csv-Mendelianr-set.csv"
cap_verb="Data from PubMed search for Mendelian randomi[s/z]ation, title only. "
xlab="Date"
ylab_root="New PubMed Entries"
mycol="dodgerblue3"
lwidth=.4
outpath='data'
outfroot="mrlit_pubmed"
outfile=paste0(outfroot,"_",gsub("-","",Sys.Date()),".csv")

## Load data ####
rdat <- read_csv(file.path(fpath, fname))
max_pubmed_date=max(rdat$`Create Date`,na.rm=T)

time_levels=c("weeks","months","years")


## Time-wise summary ####

## loop through time levels
res=bind_rows(lapply(c(time_levels),function(z){
  
  ## assign time groups then remove record duplications from e.g. pub comments
  dat=rdat|>
    mutate(tend=ceiling_date(`Create Date`,unit=z)-1)|> 
    group_by(Title,Authors)|>
    filter(`Create Date`==min(`Create Date`))|>
    ungroup()
  
  ## ensure all date rows exist even as zero.
  agg=left_join(crossing(
    tend=seq.Date(from=min(dat$tend),to=max(dat$tend),by=paste(1,gsub("s$","",z))))|>
      mutate(tend=ceiling_date(tend,unit=z)-1),
    dat|>
      group_by(tend)|>
      summarise(n=as.numeric(n()),.groups='drop'),by=join_by(tend))|>
    mutate(n=if_else(!is.finite(n),as.numeric(0),n))|>
    ungroup()|> 
    mutate(time_level=z) 
  agg
}))|> 
  select(pubmed_date=tend,n_publications=n,time_level) |> 
  mutate(cap_verb=cap_verb,file_name=fname) |> 
  mutate(time_level=stringr::str_to_title(time_level)) |> 
  mutate(last_pubmed_date=max_pubmed_date)

## write ####
# write_csv(res,paste0(outpath,"/",outfile))


## plot ####
p=ggplot(res|> 
           filter(time_level == "Weeks" & pubmed_date!=max(pubmed_date,na.rm=F)),aes(pubmed_date,n_publications))+
  # geom_smooth(se=F,colour='black')+
  geom_point(size=.25,colour=mycol,alpha=1)+
  geom_line(alpha=.75,colour=mycol)+
  geom_line(data=res %>% filter(time_level == "Weeks"),linetype=2,alpha=.75,colour=mycol)+
  scale_x_date(date_minor_breaks="1 year") +
  
  labs(x="Date",
       y=paste(ylab_root,"per week"),
       caption=paste("Most recent publication date:", max(subset(res, time_level=="Weeks")$pubmed_date)),
       x="Date") +
  theme_mrlit()
p
d <- max(subset(res, time_level=="Weeks")$pubmed_date)

ggsave(p, file=paste0("mr-weeks-", d, ".pdf"), width=8, height=4.5)
