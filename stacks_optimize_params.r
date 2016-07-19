
#input files 
#each vcf is in a folder called stacks<#> and named batch_<#>.recode.vcf)
dirpath="/borevitz/megan/PATH/TO/DIR"
vcffile_nums=c(1:64)

#count header lines in command line (grep -c "##" batch_1.vcf)
skipr=8
#non sample columns
skipc=9

#replicate sample ID numbers
rep1=c(45133:45139)
rep2=c(45373:45379)

#set output matrix
final_av=matrix(ncol=6,nrow=0)
colnames(final_av)<-c("batch","num_loci","NAs","geno_match","geno_mismatch","score")

#loop over the vcf files
for(i in vcffile_nums)
    {
        #create matrix for file data
        rep_diff=matrix(ncol=4,nrow=length(rep1))
    
        #read in vcf, skipping ##header rows
        vcf=read.delim(paste(dirpath,"/stacks",i,"/batch_",i,".recode.vcf",sep=""), sep="\t", header=T, skip=8)
        #remove non sample columns
        vcf=vcf[c(-skipc:-1)]
        #get just genotypes
        genos=apply(vcf,2,substr,1,3)
        #transpose 
        genos=t(genos)
        #put in matrix format
        genos[genos == "./."] <- "NA"
        genos[genos == "1/1"] <- "0"
        genos[genos == "0/1"] <- "1"
        genos[genos == "1/0"] <- "1"
        genos[genos == "0/0"] <- "2"
        genos.filt=genos

        #if loci and individuals remain
        if(dim(genos.filt)[1]!=0 && dim(genos.filt)[2]!=0)
            {
                #loop over each replicate
                for(r in 1:length(rep1))
                    {
                        #get data for the pair and compare
                        rep1genos=subset(genos.filt, row.names(genos.filt)==paste("trim_qc_S",rep1[r],"_il",sep=""))
                        rep2genos=subset(genos.filt, row.names(genos.filt)==paste("trim_qc_S",rep2[r],"_il",sep=""))
                        comps=suppressWarnings(as.numeric(rep1genos)==as.numeric(rep2genos))
                        
                        #calculate values
                        if(length(comps)>0)
                            {
                                rep_diff[r,1]=length(comps) #number of loci
                                rep_diff[r,2]=sum(is.na(comps)) #NAs in either or both
                                rep_diff[r,3]=sum(comps,na.rm=T) #same genos
                                rep_diff[r,4]=rep_diff[r,1]-rep_diff[r,2]-rep_diff[r,3] #diff genos
                            }else {rep_diff[r,1:4]=NA}
                    }
            }
    
        #average and add to summary
        final_av=rbind(final_av,cbind(as.integer(i), 
                        as.integer(mean(rep_diff[,1],na.rm=T)), 
                        mean(rep_diff[,2],na.rm=T), 
                        mean(rep_diff[,3],na.rm=T), 
                        mean(rep_diff[,4],na.rm=T),
                        mean(rep_diff[,3],na.rm=T)/(mean(rep_diff[,3],na.rm=T)+mean(rep_diff[,4],na.rm=T))-
                            mean(rep_diff[,2],na.rm=T)/mean(rep_diff[,1],na.rm=T)))
    }
final_av
write(t(final_av),file="stacks_params.output")

max(as.numeric(final_av[,6]), na.rm=T)
best=which(final_av[,6]==max(as.numeric(final_av[,6]), na.rm=T), arr.ind=T)
print(paste("optimal batch is",final_av[best,1]))
