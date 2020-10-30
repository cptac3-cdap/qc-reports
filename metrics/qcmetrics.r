
		options(warn = -1)
		suppressMessages(library(MSnbase))

		# figure out the directory containing this script
		# when run by Rscript

		cmdArgs <- commandArgs(trailingOnly=FALSE)
		needle <- "--file="
		match <- grep(needle, cmdArgs)
		this.dir <- dirname(sub(needle, "", cmdArgs[match]))

		# library(MSnbase)
		args <- commandArgs(TRUE)

		## this version of code accepts one file name for mzML (.gz format) and the second file name for mzID (.txt format by grep) as arguments containing all the fractions of mzML files

		mzml.file=args[1] 
		mzid.file=args[2] 
		out.file=args[3]
	
			mzml.data <- openMSfile(mzml.file) #gzfile(mzml.file)
			ms1.ms2.index.matrix.2=matrix(0, length(mzml.data),2)
			k=1
			ms1.level=-1
			for(n in 1:length(mzml.data)){
				spectrum.header=header(mzml.data,n)
				ms.level=spectrum.header$msLevel
				seqNum=spectrum.header$seqNum
				precursorScanNum=spectrum.header$precursorScanNum
				if(ms.level!=1){
					ms1.ms2.index.matrix.2[k,]=c(precursorScanNum, seqNum)
					k=k+1
				}
			}
			ms1.ms2.index.matrix.2=ms1.ms2.index.matrix.2[1:(k-1),]
	
			# Get the identified scans
			source(file.path(this.dir,"parsemzid.R"))

                        mzid <- parsemzid(mzid.file)

			# CPTAC PSMs have already been filtered at
			# 1% spectral FDR and all but the rank 1 PSMs
			# removed. In other scenarios, this filtering
			# may need to be applied

			# a scan may be mentioned more than once due to ties...
			scan.ids <- unique(mzid$psms$scan)

			# extract the rank 1 peptide sequences...
			# where there are ties we don't want to count as extra sequences
			# perhaps we should normalize the psms to sort by
			# sequence so the first PSM per scan is consistent
			# when there are ties..

			peptide.sequences <- unique(mzid$psms[mzid$psms$ordinal == 1,]$sequence)
			peptide.count <- length(peptide.sequences)
			peptide.redundancy <- length(scan.ids)/peptide.count
			f <- data.frame(table(mzid$psms[mzid$psms$ordinal == 1,]$sequence))
            if (length(rownames(f)) > 0) {
			peptide.frequency <- f[f$Freq>0,]$Freq
            } else {
            peptide.frequency <- c(0)
            }
			
			## m## none of the ms1 in mzML is found in mzID
			# ms2.per.ms1.count=c()
			# for(ms1.id in unique(ms1.ms2.index.matrix.2[,1])){
			#	ms2.ids=ms1.ms2.index.matrix.2[which(ms1.ms2.index.matrix.2[,1]==ms1.id),2]
			#	ms2.mzid.count=length(which(unique(scan.ids)%in%ms2.ids))
			#	if(ms2.mzid.count>0){
			#		ms2.per.ms1.count=c(ms2.per.ms1.count, ms2.mzid.count)
			#	}
			# }
			# ms2.per.ms1=ms2.per.ms1.count			
			
			mzml.table <- matrix(0,length(mzml.data),length(header(mzml.data))) #build a mzml.data-sized table full of zero
			#as.numeric(unlist(header(mzml.data,1))) #get the numeric value of the first row
			#names(unlist(header(mzml.data,1))) #get the column names
			for(i in 1:length(mzml.data)){
			  row.vals = as.numeric(unlist(header(mzml.data,i)))
			  mzml.table[i,] = row.vals
			} #extract the value of each row from mzml.data and put them into mzml.table(the blank matrix)
			mzml.table <- as.data.frame(mzml.table) #convert matrix to data.frame to add column & row name
			# head(mzml.table) #get the first six rows of the data frame(table)
			colnames(mzml.table) <- names(unlist(header(mzml.data,1))) #name the columns
			
			
			ms1.all.ids=as.vector(mzml.table[which(mzml.table$msLevel==1), 'seqNum'])
			ms2.per.ms1.both.table=matrix(0, length(ms1.all.ids), 3)
			k=1
			for(ms1 in ms1.all.ids){
			  mzml.ms2=as.vector(mzml.table[which(mzml.table$precursorScanNum==ms1), 'seqNum'])
			  if(length(mzml.ms2)==0){
			    ms2.per.ms1.both.table[k, ]=c(ms1, 0 ,0)
			  }else{
			    ms2.mzid.count=which(unique(scan.ids)%in%mzml.ms2)
			    ms2.per.ms1.both.table[k, ]=c(ms1, length(mzml.ms2) , length(ms2.mzid.count))
			  }
			  k=k+1
			}
			ms2.per.ms1.count=as.numeric(ms2.per.ms1.both.table[,3])
			ms2.per.ms1=ms2.per.ms1.count
			
			
			ms2.level.equal.max=length(which(ms2.per.ms1.count==max(ms2.per.ms1.count)))/length(ms2.per.ms1.count)
			
			run=runInfo(mzml.data)
			spectra.count=length(mzml.data) ## will be replaced

                        mzml.file.basename=basename(mzml.file)

                        # Assumes there is only one period in the filename...
                        spectrum.file.basename=strsplit(mzml.file.basename, ".", fixed=TRUE)[[1]][1]
                        base.name.strsplit=strsplit(spectrum.file.basename, "_")[[1]]

                        # Does this handle double underscore correctly? Hopefully...
                        analytical.sample.name=paste(base.name.strsplit[1:(length(base.name.strsplit)-1)], collapse="_")
                        # Picking out the fraction number is hard because
                        # there not much in the way of convention here,
                        # there are a bunch of possibilities. 
                        #
                        # Assume it is the number (or A) looking thing
                        # with or without a non-numeric prefix in the last
                        # "_" delimited chunk. Add other special cases as needed...

                        last_chunk=base.name.strsplit[length(base.name.strsplit)]
                        fraction.num=type.convert(gsub("^[^0-9]*([0-9]+|[A-Z])$","\\1",last_chunk),as.is=TRUE)

			## mzml.data -> unique(scan.ids)
			info.lists=c("msLevel", "precursorCharge", "precursorIntensity", "precursorMZ","retentionTime")
			info.matrix=matrix(0, length(unique(scan.ids)), length(info.lists))
			for(i in 1:length(info.lists)){
				info.list=info.lists[i]
				info.matrix[,i]=unlist(lapply(unique(scan.ids), function(n){spectrum.header=header(mzml.data,n);as.numeric(spectrum.header[info.list])}))
			}
			precursorMW=(info.matrix[,4]-1.0078)*info.matrix[,2]
			spectra.count=length(unique(scan.ids))

			result.info=c(
				spectrum.file.basename,
				analytical.sample.name,
				fraction.num,
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 3], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 3]),
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 4], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 4]),				
				as.numeric(quantile(precursorMW[which(info.matrix[,1]==2)], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(precursorMW[which(info.matrix[,1]==2)]),				
				quantile(ms2.per.ms1, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)),
				mean(ms2.per.ms1),				
				length(info.matrix[which(info.matrix[,1]==2),1]),
                                length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==0), 1]),
                                length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==1), 1]),
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==2), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==3), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==4), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]>4), 1]), 
				ms2.level.equal.max,
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 5], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 5])
			)
			result.info=data.frame(t(result.info))
			colnames(result.info)=c("spectrumBasename", "analyticalBasename", "fractionNum", 
									"minPrecursorIntensity", "PrecursorIntensity5perc", 
									"PrecursorIntensity25perc", "PrecursorIntensity50perc", "PrecursorIntensity75perc", 
									"PrecursorIntensity95perc", "maxPrecursorIntensity", "meanPrecursorIntensity",
									"minPrecursorMZ", "PrecursorMZ5perc", "PrecursorMZ25perc", "PrecursorMZ50perc", 
									"PrecursorMZ75perc", "PrecursorMZ95perc", "maxPrecursorMZ", "meanPrecursorMZ",
									"minPrecursorMW", "PrecursorMW5perc", "PrecursorMW25perc", "PrecursorMW50perc", 
									"PrecursorMW75perc", "PrecursorMW95perc", "maxPrecursorMW", "meanPrecursorMW",
									"minMS2perMS1", "MS2perMS1_5perc", "MS2perMS1_25perc", "MS2perMS1_50perc", 
									"MS2perMS1_75perc", "MS2perMS1_95perc", 
									"maxMS2perMS1", "meanMS2perMS1",
									"numofMS2", 
									"MS2Charge0", "MS2Charge1", "MS2Charge2", "MS2Charge3", "MS2Charge4", "MS2Charge5andHigher", "ratioMaxMS2Count",
									"minRetentionTime", "RetentionTime5perc", 
									"RetentionTime25perc", "RetentionTime50perc", "RetentionTime75perc", 
									"RetentionTime95perc", "maxRetentionTime", "meanRetentionTime"
									)
									
		mzml.id= result.info
	

			
		############## part 2  on mzml data only


			
		## this version of code accepts one file name as argument in .gz (or uncompressed) format
		
			#mzml.data <- openMSfile(mzml.file)			
			
		  run=runInfo(mzml.data)
			spectra.count=length(mzml.data)
			mzml.file.basename=basename(mzml.file)
			
			# Assumes there is only one period in the filename...
			spectrum.file.basename=strsplit(mzml.file.basename, ".", fixed=TRUE)[[1]][1]
			base.name.strsplit=strsplit(spectrum.file.basename, "_")[[1]]
			
			# Does this handle double underscore correctly? Hopefully...
			analytical.sample.name=paste(base.name.strsplit[1:(length(base.name.strsplit)-1)], collapse="_")
			# Picking out the fraction number is hard because
			# there not much in the way of convention here,
			# there are a bunch of possibilities. 
			#
			# Assume it is the number (or A) looking thing
			# with or without a non-numeric prefix in the last
                        # "_" delimited chunk. Add other special cases as needed...

			last_chunk=base.name.strsplit[length(base.name.strsplit)]
			fraction.num=type.convert(gsub("^[^0-9]*([0-9]+|A)$","\\1",last_chunk),as.is=TRUE)
			
			info.lists=c("msLevel", "precursorCharge", "precursorIntensity", "precursorMZ","retentionTime")
			info.matrix=matrix(0, length(mzml.data), length(info.lists))
			for(i in 1:length(info.lists)){
				info.list=info.lists[i]
				info.matrix[,i]=unlist(lapply(1:length(mzml.data), function(n){spectrum.header=header(mzml.data,n);as.numeric(spectrum.header[info.list])}))
			}
			
			# comment out the old way of getting MS2 per MS1 for mzML file
			#level.string=paste(info.matrix[,1], collapse="")
			#level.string.vector=unlist(strsplit(level.string, "1"))
			#level.string.count=unlist(lapply(level.string.vector, function(v){nchar(v)}))
			
			# get the data directly from part I, need the same mzML file
			level.string.count=as.numeric(ms2.per.ms1.both.table[,2])
			# level.string.count=level.string.count[which(level.string.count>0)] # added 2/23/2018 # then delete 2/28/2018
			
			ms2.level.equal.max=length(which(level.string.count==max(level.string.count)))/length(level.string.count)
			precursorMW=(info.matrix[,4]-1.0078)*info.matrix[,2]
			
			
			ms2.per.ms1=level.string.count	
			result.info=c(
				spectrum.file.basename,
				analytical.sample.name,
				fraction.num,
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 3], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 3]),
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 4], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 4]), 
				as.numeric(quantile(precursorMW[which(info.matrix[,1]==2)], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(precursorMW[which(info.matrix[,1]==2)]),				
				quantile(ms2.per.ms1, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)),
				mean(ms2.per.ms1),			
				length(info.matrix[which(info.matrix[,1]==2),1]),
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==0), 1]),
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==1), 1]),
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==2), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==3), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]==4), 1]), 
				length(info.matrix[which(info.matrix[,1]==2 & info.matrix[,2]>4), 1]), 
				ms2.level.equal.max,
				as.numeric(quantile(info.matrix[which(info.matrix[,1]==2), 5], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))),
				mean(info.matrix[which(info.matrix[,1]==2), 5])
			)
			result.info=data.frame(t(result.info))
			colnames(result.info)=c("spectrumBasename", 
				"analyticalBasename", "fractionNum", 
													
				"minPrecursorIntensity", "PrecursorIntensity5perc", 
				"PrecursorIntensity25perc", "PrecursorIntensity50perc", 
				"PrecursorIntensity75perc", 
													
				"PrecursorIntensity95perc", "maxPrecursorIntensity", 
				"meanPrecursorIntensity",
													
				"minPrecursorMZ", "PrecursorMZ5perc", "PrecursorMZ25perc", 
				"PrecursorMZ50perc", "PrecursorMZ75perc", "PrecursorMZ95perc", 
				"maxPrecursorMZ", "meanPrecursorMZ",
													
				"minPrecursorMW", "PrecursorMW5perc", "PrecursorMW25perc", 
				"PrecursorMW50perc", "PrecursorMW75perc", "PrecursorMW95perc", 
				"maxPrecursorMW", "meanPrecursorMW",
													
				"minMS2perMS1", "MS2perMS1_5perc", "MS2perMS1_25perc", 
				"MS2perMS1_50perc", "MS2perMS1_75perc", "MS2perMS1_95perc", 
				"maxMS2perMS1", "meanMS2perMS1",
													
				"numofMS2", 
													
				"MS2Charge0", "MS2Charge1", "MS2Charge2", "MS2Charge3", "MS2Charge4", 
				"MS2Charge5", "ratioMaxMS2Count",
									"minRetentionTime", "RetentionTime5perc", 
									"RetentionTime25perc", "RetentionTime50perc", "RetentionTime75perc", 
									"RetentionTime95perc", "maxRetentionTime", "meanRetentionTime"
									
			)

		mzml.all= result.info



	mzml.all.id=matrix(0, dim(mzml.id)[1], 3+(ncol(mzml.id)-3)*2)
	colnames(mzml.all.id)=1:dim(mzml.all.id)[2]
	mzml.all.id[,1]=as.vector(mzml.all[,1])
	mzml.all.id[,2]=as.vector(mzml.all[,2])
	mzml.all.id[,3]=as.vector(mzml.all[,3])
	k=4
	colnames(mzml.all.id)[1:3]=colnames(mzml.all)[1:3]
	for(i in 4:ncol(mzml.all)){
		mzml.all.id[, k]=as.vector(mzml.all[, i])
		mzml.all.id[, k+1]=as.vector(mzml.id[, i])
		colnames(mzml.all.id)[k]=paste(colnames(mzml.all)[i], "ALL", sep="_")
		colnames(mzml.all.id)[k+1]=paste(colnames(mzml.all)[i], "ID", sep="_")
		k=k+2
	}

	extra = data.frame(PeptideCount=c(peptide.count),
			   PeptideRedundancy=c(peptide.redundancy))
	mzml.all.id <- cbind(mzml.all.id,extra)

	pepfreq <- as.numeric(quantile(peptide.frequency, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))
	pepfreq <- data.frame(matrix(c(pepfreq, mean(peptide.frequency)),nrow=1))
	colnames(pepfreq) <- c("minPeptideRedundancy", "PeptideRedundancy5perc",
                               "PeptideRedundancy25perc", "PeptideRedundancy50perc", "PeptideRedundancy75perc",
                               "PeptideRedundancy95perc", "maxPeptideRedundancy", "meanPeptideRedundancy")
	mzml.all.id <- cbind(mzml.all.id,pepfreq)
	mzml.all.id[is.na(mzml.all.id)] <- 0
	write.table(mzml.all.id, out.file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
